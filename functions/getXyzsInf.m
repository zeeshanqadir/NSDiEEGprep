%% This function returns inflated gifti XYZ positions from an electrodes.tsv table
% Adapted from DH ieeg_snap2inflate.m to perform separately for each hemisphere
%
%   xyzsInf = getXyzsInf(electrodes, hemi, gii, giiInf);
%   xyzsInf = getXyzsInf(electrodes, hemi, gii, giiInf, distThresh);
%       electrodes =        nx_ table, read from an electrodes.tsv file. Must contain column "hemisphere" (str, 'R' or 'L').
%                               Must also contain columns "Destrieux_label" (numerical) and "Destrieux_label_text" if <checks> is given as true.
%       hemi =              char, 'r' or 'l', corresponding to the hemisphere of the pial and inflated giftis
%       gii =               gifti object, pial surface of one hemisphere
%       giiInf =            gifti object, inflated surface of one hemisphere
%       distThresh =        double (optional), Only electrodes within <distThresh> mm of a pial vertex are kept in output.
%                               Default = 5.
%       checks =            logical (optional), whether to more strictly check that each electrode is cortical and not WM/thalamic/outside brain. Default == false
%
%   Returns:
%       xyzsInf =           nx3 double, XYZ positions of gray matter electrodes in inflated gifti space. Electrodes in
%                               white matter, contralateral hemisphere, or beyond <distThresh> are returned as [NaN, NaN, NaN] rows.
%
%   HH 2021
%
function xyzsInf = getXyzsInf(electrodes, hemi, gii, giiInf, distThresh, checks)

    if nargin < 6, checks = false; end % to perform checks based on Destrieux atlas names
    if nargin < 5 || isempty(distThresh), distThresh = 5; end % maximum distance allowed between valid electrode and gifti vertex
    
    hemi = lower(hemi);
    assert(ismember(hemi, {'r', 'l'}), "hemi must be given as 'r' or 'l'");
    assert(size(gii.vertices, 1) == size(giiInf.vertices, 1), 'gii and giiInf do not have the same number of vertices');
    
    xyzs = [electrodes.x, electrodes.y, electrodes.z];
    hemiLabs = lower(electrodes.hemisphere);

    if checks
        anatLabs = electrodes.Destrieux_label; % numerical labels
        anatText = electrodes.Destrieux_label_text; % char labels
    else
        anatLabs = zeros(height(electrodes), 1);
        anatText = cell(height(electrodes), 1);
    end

    xyzsInf = nan(size(xyzs));
    for ii = 1:size(xyzs, 1)
        
        if ~strcmp(hemiLabs(ii), hemi) % wrong hemisphere, continue
            continue
        end
        
        if checks && (isnan(anatLabs(ii)) || ismember(anatLabs(ii), [0, 41, 49, 10])) % non existent, outside brain, WM, or L/R thalamic
            continue
        end
        
        if checks && (~startsWith(anatText{ii}, sprintf('%sh', hemi))) % must start with lh or rh to be a cortical site
            continue
        end
        
        dists = vecnorm(gii.vertices - xyzs(ii, :), 2, 2); % distance from each vertex in gifti to current electrode
        [minDist, idx] = min(dists);
        
        if minDist <= distThresh
            xyzsInf(ii, :) = giiInf.vertices(idx, :);
        end
        
    end
    
end