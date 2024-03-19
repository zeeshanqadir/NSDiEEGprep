%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes SOZ

%% Paths, load tables

% cd here
setPathsNSD;
addpath('functions');

sub = 'X';
[tasks, runs] = subjectConfig(sub);

% electrodes path (for SOZ)
elecPath = fullfile(localRoot, sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub));
elecs = ieeg_readtableRmHyphens(elecPath);

% load file to indicate which tasks were performed in which sessions (to use as a possible random effect later)
sessions = readtable(fullfile('data', 'metadata', sprintf('%s_sessions.txt', sub)), 'FileType', 'text', 'Delimiter', '\t');

% All possible runs: col1 = mefPath, col2 = channels path, col3 = events path
dataPathStems = {};
taskruns = [];
for ii = 1:length(runs)
    for jj = 1:length(tasks)
        dataPathStems = [dataPathStems;
            {fullfile(localRoot, sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_task-NSDspecial%02d_run-%02d', sub, tasks(jj), runs(ii)))}];
            taskruns = [taskruns; tasks(jj), runs(ii)]; % keep track of which NSD task and run
    end
end

% convert taskruns to table, this way 3rd column later will be for bad channels
taskruns = array2table(taskruns, 'VariableNames', {'tasknumber', 'run'});
taskruns.badChs = repmat({'n/a'}, height(taskruns), 1);

% contants
srate = 1200;

%% Loop through NSD01 - NSD10 and all runs individually, concatenate output at the end

Mephys = []; % ephys data matrix
Meye = []; % eyetracker data matrix
eventsST = []; % stim/test events

for ii = 1:length(dataPathStems)
    %% Load channels, events, mef data

    fprintf('Loading %s\n', dataPathStems{ii});

    % load channels and events
    try
        channels = ieeg_readtableRmHyphens(sprintf('%s_channels.tsv', dataPathStems{ii}));
        events =  readtable(sprintf('%s_events.tsv', dataPathStems{ii}), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');
    catch
        warning('Error loading channels/events file for %02d, skipping', ii); % e.g., no 2nd run for a certain NSD session
        continue
    end

    % force status_description to be text (so that it can be searched) in case all nans
    events.status_description = string(events.status_description);

    % load mef data
    [metadata, data] = readMef3(sprintf('%s_ieeg.mefd', dataPathStems{ii}));

    %% Reorganize and clean up events so that adjacent rest and stim events are combined

    % First event should be rest. Sometimes this pdiode onset is not acquired bc pdiode was on during presequence gray
    if strcmp(events.trial_type(1), 'stim')
        warning('Appending artificial rest trial to start of events');
        T1 = cell2table({events.onset(1) - 2, 2, events.sample_start(1) - 2*srate, 'rest', '', 'good', 'n/a'}, ...
                        'VariableNames', events.Properties.VariableNames);
        events = [T1; events];
    end

    % reassemble events table to be stim/test only. Add pre (rest) status and durations
    idxST = find(ismember(events.trial_type, {'stim', 'test'}));
    assert(all(strcmp(events.trial_type(idxST-1), 'rest')), 'Error: trials preceding some stim/test trials are not rest');
    eventsSTCurr = events(idxST, :);
    eventsSTCurr.pre_duration = events.duration(idxST-1, :);
    eventsSTCurr.pre_status = events.status(idxST-1, :);
    eventsSTCurr.pre_status_description = events.status_description(idxST-1, :);

    % delete events that aren't good currently or in pre-rest, as well as high interictal
    eventsSTCurr(~strcmp(eventsSTCurr.status, 'good') | ~strcmp(eventsSTCurr.pre_status, 'good'), :) = [];
    eventsSTCurr(contains(eventsSTCurr.status_description, 'Categorical-level/High') | contains(eventsSTCurr.pre_status_description, 'Categorical-level/High'), :) = [];

    % add random effects time markers: task (NSD1-10), run (1, 2...), session (ordered blocks of time)
    eventsSTCurr.tasknumber = repmat(taskruns.tasknumber(ii), height(eventsSTCurr), 1);
    eventsSTCurr.run = repmat(taskruns.run(ii), height(eventsSTCurr), 1);
    sess = sessions.session(sessions.task == taskruns.tasknumber(ii) & sessions.run == taskruns.run(ii));
    eventsSTCurr.session = repmat(sess, height(eventsSTCurr), 1);

    %% Preprocess and convert to trial struction
    
    % Highpass the SEEG channels
    data(strcmp(channels.type, 'SEEG'), :) = ieeg_highpass(data(strcmp(channels.type, 'SEEG'), :)', srate)';

    % Convert data to trial structure
    samprange = -2*srate:2*srate-1;

    % remove any events missing the full sample range (e.g. due to crash)
    eventsSTCurr(eventsSTCurr.sample_start+samprange(end)+1 > size(data, 2), :) = [];

    tt = samprange/srate;
    M = zeros(size(data, 1), length(samprange), height(eventsSTCurr));
    for jj = 1:height(eventsSTCurr)
        M(:, :, jj) = data(:, (eventsSTCurr.sample_start(jj)+samprange(1)+1) : (eventsSTCurr.sample_start(jj)+samprange(end)+1)); % convert to 1 indexing
    end

    % check the DI1 channel to make sure that the first sample of 1 occurs at t=0
    %di1 = squeeze(M(strcmp(channels.name, 'DigitalInput1'), :, :));
    %figure; plot(tt, di1); xlim([-0.01, 0.01]);

    % Separate ephys from eyetracker data
    MephysCurr = M(strcmp(channels.type, 'SEEG'), :, :);
    if ii > 1
        assert(sum(strcmp(channels.type, 'SEEG')) == height(channelsEphys), 'Error: Mismatch between current number of ephys channels and prev number of ephys channels');
    end
    channelsEphys = channels(strcmp(channels.type, 'SEEG'), :);
    MeyeCurr = M(startsWith(channels.name, 'Eyetracker'), :, :); assert(size(MeyeCurr, 1) == 14, 'Error: Should have 14 eyetracker channels');
    eyetrackerLabels = channels.name(startsWith(channels.name, 'Eyetracker'));

    % Find all bad channels (including SOZ)
    elecsSOZ = elecs.name(contains(elecs.seizure_zone, 'SOZ')); % first set SOZ channels to bad
    channelsEphys.status(ismember(channelsEphys.name, elecsSOZ)) = {'bad'};
    badChs = union(ieeg_getChsNR(channelsEphys, elecs), find(~strcmp(channelsEphys.status, 'good'))); % double check to see if any NR channels have been missed
    
    % CAR by lowest variance, with all trials considered the same group
    opts = struct();
    opts.vartype = 'var';
    [MephysCurr, chsUsed] = ccep_CARVariance(tt, MephysCurr, srate, badChs, [], opts);

    %% Concatenate data, add info on good channels at current run

    Mephys = cat(3, Mephys, MephysCurr);
    Meye = cat(3, Meye, MeyeCurr);
    eventsST = [eventsST; eventsSTCurr];
    taskruns.badChs{ii} = badChs;

end

% remove taskruns entries that didn't happen
taskruns(strcmp(taskruns.badChs, 'n/a'), :) = [];

%% Save all outputs

% use most recent channelsephys (all are same height bc assertion), but remove status column. bad channels are indicated in taskruns variable
if any(ismember(channelsEphys.Properties.VariableNames, {'status', 'status_description'}))
    channelsEphys.status = [];
    channelsEphys.status_description = [];
end

outdir = fullfile('output', 'fromMef', sub);
mkdir(outdir);
save(fullfile(outdir, sprintf('preprocessed_%s.mat', sub)), 'tt', 'srate', 'Mephys', 'Meye', 'eventsST', 'channelsEphys', 'taskruns', 'eyetrackerLabels', '-v7.3');

writetable(eventsST, fullfile(outdir, sprintf('eventsSTPreproc_%s.tsv', sub)), 'FileType', 'text', 'Delimiter', '\t');

% Plot and save a heatmap of ERPs, for channels good in at least 5 runs
% nan-out trials for channels in bad runs before averaging and plotting
MephysPlot = Mephys;
for ii = 1:height(taskruns)
    MephysPlot(taskruns.badChs{ii}, :, eventsST.tasknumber == taskruns.tasknumber(ii) & eventsST.run == taskruns.run(ii)) = nan;
end
evoked = mean(MephysPlot, 3, 'omitnan'); % average across trials for each channel (only for runs where channels are good)
[goodChsOnN, badChsOnN] = goodChsNRuns(taskruns.badChs, height(channelsEphys), 5); % must be good in 5 runs to be a good channel to plot
evoked(badChsOnN, :) = nan;
figure('Position', [200, 200, 600, 1200]);
imagescNaN(tt, 1:height(channelsEphys), evoked, 'CLim', [-50, 50]);
xline(0, 'Color', 'r');
xlim([-0.8, 0.8]);
set(gca, 'YTick', 1:2:height(channelsEphys), 'YTickLabels', channelsEphys.name(1:2:end));
title('ERPs');
xlabel('Time (s)'); ylabel('Channels');
colorbar;
saveas(gcf, fullfile(outdir, sprintf('ERPsummary_%s.png', sub)));
