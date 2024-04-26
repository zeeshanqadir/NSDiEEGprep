%% This script calculates MNI152 and MNI305 positions each subject and saves them

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

ss = 17;
sub_label = sprintf('%02d', ss);

ses_label = 'ieeg01';
outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);

% load electrodes
elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);
elecmatrix = [elecs.x, elecs.y, elecs.z];

%% Get and save MNI152 positions to electrodesMni152.tsv (volumetric, SPM12)

% locate forward deformation field from SPM. There are variabilities in session name, so we use dir to find a matching one
niiPath = dir(fullfile(localDataPath.input, 'sourcedata', 'spm_forward_deformation_fields', sprintf('sub-%s_ses-*_T1w_acpc.nii', sub_label)));
assert(length(niiPath) == 1, 'Error: did not find exactly one match in sourcedata T1w MRI'); % check for only one unique match
niiPath = fullfile(niiPath.folder, niiPath.name);

% create a location in derivatives to save the transformed electrode images
rootdirMni = fullfile(localDataPath.input, 'derivatives', 'MNI152_electrode_transformations', sprintf('sub-%s', sub_label));
mkdir(rootdirMni);

% calculate MNI152 coordinates for electrodes
xyzMni152 = ieeg_getXyzMni(elecmatrix, niiPath, rootdirMni);

% save as separate MNI 152 electrodes table
elecsMni152Path = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_space-' 'MNI152NLin2009' '_electrodes.tsv']);
elecsMni152 = elecs;
elecsMni152.x = xyzMni152(:, 1); elecsMni152.y = xyzMni152(:, 2); elecsMni152.z = xyzMni152(:, 3);
writetable(elecsMni152, elecsMni152Path, 'FileType', 'text', 'Delimiter', '\t');

fprintf('Saved to %s\n', elecsMni152Path);

%% Get and save MNI305 positions (through fsaverage)

% FS dir of current subject
FSdir = fullfile(localDataPath.input, 'sourcedata', 'freesurfer', sprintf('sub-%s', sub_label));
FSsubjectsdir = fullfile(FSdir, '..');

% calculate MNI305 coordinates for electrodes
[xyzMni305, vertIdxFsavg, minDists, surfUsed] = ieeg_mni305ThroughFsSphere(elecmatrix, elecs.hemisphere, FSdir, FSsubjectsdir, 'closest', 5);

% save as separate MNI 305 electrodes table
elecsMni305Path = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_space-' 'MNI305' '_electrodes.tsv']);
elecsMni305 = elecs;
elecsMni305.x = xyzMni305(:, 1); elecsMni305.y = xyzMni305(:, 2); elecsMni305.z = xyzMni305(:, 3);
elecsMni305.vertex_fsaverage = vertIdxFsavg; % also add a column to indicate vertex on fsavg, so we can easily get position for inflated brain
writetable(elecsMni305, elecsMni305Path, 'FileType', 'text', 'Delimiter', '\t');

fprintf('Saved to %s\n', elecsMni305Path);

%% Normalize bb power per run

% Initialize normalized log power of BB
Mbb_norm = log10(Mbb); 

% Indicate the interval for baseline, used in normalization
norm_int = find(tt>-.2 & tt<0);

% Normalize per run
for run_idx = 1:max(eventsST.tasknumber)
    this_run = find(eventsST.tasknumber==run_idx); % out of 1500
    
    % find pre-stim events with 'good' status
    trials_norm = find(ismember(eventsST.pre_status,'good') & eventsST.tasknumber==run_idx);

    Mbb_norm(:,:,this_run) = minus(Mbb_norm(:,:,this_run),mean(Mbb_norm(:,norm_int,trials_norm),[2 3],'omitnan'));
end

%% Find repeated images, calculate SNR

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);

all_chan_snr = NaN(size(Mbb_norm,1),1);
t_avg = tt>0.1 & tt<.5;

for el_nr = 1:size(Mbb_norm,1)
    
    if ismember(all_channels.type(el_nr),'SEEG') && all_channels.status(el_nr)==1
        bb_strength = squeeze(mean(Mbb_norm(el_nr,t_avg==1,:),2));
        
        all_repeats = find(nsd_repeats>0);
        shared_idx_repeats = unique(shared_idx(all_repeats)); % 100 images
        repeats_bb_strength = cell(length(shared_idx_repeats),1);
        for kk = 1:length(shared_idx_repeats)
            these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
            repeats_bb_strength{kk} = bb_strength(these_trials); 
        end
        
        [NCSNR, p, NCSNRNull] = estimateNCSNR(repeats_bb_strength, 1000);
        all_chan_snr(el_nr) = NCSNR;
    end
end

%% render and plot noise ceiling SNR

elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);

name = all_channels.name;
all_channels_table = table(name);
elecs = ieeg_sortElectrodes(elecs, all_channels_table, 0);

% load pial and inflated giftis
gL = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['white.L.surf.gii']));
gR = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['white.R.surf.gii']));
gL_infl = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['inflated.L.surf.gii']));
gR_infl = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['inflated.R.surf.gii']));

% snap electrodes to surface and then move to inflated
xyz_inflated = ieeg_snap2inflated(elecs,gR,gL,gR_infl,gL_infl,4);