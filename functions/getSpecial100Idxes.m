%% Finds trial indices corresponding to each special 100 image, by path name
% Returns a 100 x 2 cell array. First column = name of special100 file, second column = trial indices
%
function special100Idxes = getSpecial100Idxes(stimFiles)

    if exist('data/stims/special100Names.txt', 'file') % read stim names directly if they've been written
        stimNames = readcell('data/stims/special100Names.txt', 'Delimiter', ',');
        assert(length(stimNames) == 100, 'Error: special100Names did not contain 100 entries');
    else % get names from image set and save as output
        % Get names of stimuli
        stimPaths = glob('data/stims/special100/*.png');
        assert(length(stimPaths) == 100, 'Error: Did not find 100 special100 files');
        stimNames = extractBetween(stimPaths, 'special100/', '.png');
        writecell(stimNames, 'data/stims/special100Names.txt');
    end

    % find indices of each special100 image
    idxes = cell(100, 1);
    for ii = 1:100
        idxes{ii} = find(contains(stimFiles, stimNames{ii}));
        assert(~isempty(idxes{ii}), 'Error: No trial index found for special100 image %d\n', ii);
    end

    % Add name and indices
    special100Idxes = [stimNames, idxes];

end