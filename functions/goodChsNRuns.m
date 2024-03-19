% This function returns which channels are good across at least n runs, where the input is bad chs for each run
%
%   goodChs = goodChsNRuns(badChs, numChs, n);
%   [goodChs, badChs] = goodChsNRuns(badChs, numChs, n);
%
%   INPUTS
%       badChsRuns =        1xn numerical cell array. bad channels in each run
%       numChs =        num. Total number of channels
%       n      =        num. Minimum number of runs for a channel to be good in
%
%   OUTPUTS
%       goodChsOnN =       1xp num. Indices of good channels across >= n runs
%       badChsOnN =        1xq num. Indices of bad channels (complement of goodChs on 1:numChs)
%
function [goodChsOnN, badChsOnN] = goodChsNRuns(badChsRuns, numChs, n)

    numRuns = length(badChsRuns);
    assert(n >= 1 && n <= numRuns, 'Error: n must be specified between 1 and %d runs', numRuns);
    
    % chs x runs. 1 indicates a bad channel, 0 indicates a good channel
    badChsLogical = zeros(numChs, numRuns, 'logical');
    for ii = 1:numRuns
        badChsLogical(badChsRuns{ii}, ii) = true;
    end

    % find which rows have at least n 0s
    goodChsLogical = sum(~badChsLogical, 2) >= n;

    % convert to numerical indices, and also return complement as bad chs
    goodChsOnN = find(goodChsLogical);
    badChsOnN = find(~goodChsLogical);

end