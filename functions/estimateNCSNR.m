%% This function estimates the NCSNR of input signal strengths, as the SD of signal divided by the SD of noise
%   Derived from estimateNoiseCeiling.m
%
%   Variance due to noise is estimated from trials of the same stimulus, averaged across all stimuli. The remaining
%   variance is explained by signal, at the single trial level. If numAvg is given, will return the noise ceiling
%   corresponding to averaging numAvg trials together (noise variance is reduced numAvg-fold, assuming Gaussian).
%
%   NCSNR = estimateNoiseCeiling(strength);
%   [NCSNR, p] = estimateNoiseCeiling(strength, nperm);
%       strength =          nx1 cell, measure of scalar signal strength for each trial with n unique stimuli and variable trials per stimulus.
%                               Examples include average bandpower or a coefficient fit for the trial on some template (e.g. fmri beta)
%       nPerm =            (optional) positive integer, number of permutations to use to construct permutation null model. Default = 1000
%
%   Returns:
%       NCSNR =             double, SNR as (SD of Signal)/(SD of noise).
%
%   Based on:
%   Allen, E.J., St-Yves, G., Wu, Y. et al. A massive 7T fMRI dataset to bridge cognitive neuroscience and artificial intelligence.
%       Nat Neurosci 25, 116?126 (2022). https://doi.org/10.1038/s41593-021-00962-x
%
%   For derivations and detailed explanation, please see noise ceiling estimation in Methods section of article.
%       - Some more info is in "Functional data (NSD)/Noise ceiling" in NSD data manual, https://cvnlab.slite.com/p/channel/CPyFRAyDYpxdkPK6YbB5R1
%
%   Demo:
%       NCSNRtest = zeros(10000, 1);
%       for ii = 1:10000
%           xx = mat2cell(rand*10*randn(600, 1) + rand*10, 6*ones(100, 1));
%           NCSNRtest(ii) = estimateNCSNR(xx);
%       end
%       figure; histogram(NCSNRtest, 30);
%
% HH 2022/05
%
function [NCSNR, p, NCSNRNull] = estimateNCSNR(strength, nPerm)


    % Formatting: vertical cell array and vertical array in each cell, remove nan entries
    strength = strength(:);
    strength = cellfun(@(x) x(:), strength, 'UniformOutput', false);
    strength = cellfun(@(x) x(~isnan(x)), strength, 'UniformOutput', false);
    
    % number of entries in each cell
    cellLengths = cellfun(@length, strength);
    assert(sum(cellLengths >= 3) > 0.5*length(cellLengths), 'Error: More than half of conditions have <3 trials. Noise SD is poorly estimated');

    % Zscore strength so SD = var = 1
    meanStrength = mean(cell2mat(strength));
    SDStrength = std(cell2mat(strength), 0); % using N-1 unbiased estimator
    Zstr = cellfun(@(x) (x - meanStrength)/SDStrength, strength, 'UniformOutput', false);
    
    % Calculate NCSNR
    NCSNR = calcNCSNR(Zstr);

    if nargout < 2, return; end % no p-value requested

    % Calculate null distribution and p-value
    if nargin < 2 || isempty(nPerm), nPerm = 1000; end
    assert(nPerm > 0, 'Error: nPerm must be a positive integer');
    NCSNRNull = nan(nPerm, 1);
    ZstrAll = cell2mat(Zstr);
    n = length(ZstrAll); % total number of strength entries
    for ii = 1:nPerm
        Perm = randperm(n);
        ZstrNull = mat2cell(ZstrAll(Perm), cellLengths);
        NCSNRNull(ii) = calcNCSNR(ZstrNull);
    end
    p = sum(NCSNRNull >= NCSNR) / nPerm;
    
end

% Helper function so as to call it directly with normalized ZstrNull
function NCSNR = calcNCSNR(Zstr)

    cellLengths = cellfun(@length, Zstr);

    % residual variance across trials of each stimuli, averaged across stimuli.
    % Reason to avg stimuli with >= 3 trials bc Bessel correction is still biased for SD (underestimates noise SD here). 3 matches sample size in NSD fMRI
    varNoise = mean(cellfun(@(x) var(x, 0), Zstr(cellLengths >= 3)));
    sdNoise = sqrt(varNoise);
    
    varSig = 1 - varNoise; % variability explained by signal. Check this distribution using demo. White noise is expected to have varSig ~0 on avg
    
    if varSig <= 0
        sdSig = 0; % half-rectification of signal SD (signal variance can be negative like R^2)
    else
        sdSig = sqrt(varSig); % sqrt of variance explained by signal
    end
    NCSNR = sdSig/sdNoise; % noise ceiling signal-to-noise ratio

end