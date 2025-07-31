%% Function to plot most/least responsive NSD images

function im_idx = plotResponsiveImages( data, shared_idx, im_path, varargin)
 
% input:
%     data: iEEG/fMRI data vector (n x 1)
%     shared_idx: selected NSD shared ID from pool of 1000  (n x 1); 1 >= n >= 1000
%     im_path: path to load NSD shared 1000 images from
%     varargin{1}: 'max' for most responsive (default) or 'min' for least responsive
%     varargin{2}: No. of images to display, i (default: i=25)
%     varargin{3}: subset NSD shared ID to highlight j of i; j<=i (see use case) 
%
% output:
%     im_idx: sorted order of n images
%
% Example usage:
%
%   1. Display 200 most responsive images of shared 1000
%   tempbb  =  mean( bb_el_1000(winresp, :), 1, 'omitnan')';
%   im_path = '..data/nsddata_stimuli/stimuli/shared1000/';
%   im_idx  = plotResponsiveImages( tempbb, 1:1000, impath, 'max', 200); 
%
%   2. Display 200 most responsive images of shared 1000 highlighting all faces
%   tempbb  =  mean( bb_el_1000(winresp, :), 1, 'omitnan')';
%   im_path = '..data/nsddata_stimuli/stimuli/shared1000/';
%   im_idx  = plotResponsiveImages( tempbb, 1:1000, impath, 'max', 200, find(face_idx)); 
% 
%   3. Display 50 most responsive images of special100
%   tempbb  =  mean( bb_el_1000(winresp, find(special100_idx)), 1, 'omitnan')';
%   im_path = '..data/nsddata_stimuli/stimuli/shared1000/';
%   im_idx  = plotResponsiveImages( tempbb, find(special100_idx), impath, 'max', 50); 
%
% ZQ Nov 2024


    p = inputParser;

    addRequired(p, 'data', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'shared_idx', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'im_path', @ischar);

    parse(p, data, shared_idx, im_path);

    data = p.Results.data;
    shared_idx = p.Results.shared_idx;
    im_path = p.Results.im_path;

    if size(data,2)>1,   data= data';  end
    if size(shared_idx,2)>1,   shared_idx= shared_idx';  end

    assert(size(data,1) == size(shared_idx,1), 'Data and shared indices list should have same number of rows (samples)!')

    if isempty(im_path)
        error('Path to image directory missing!')
    end

    if nargin == 3
        
        inSt = 'max';
        inInt = min(length(data), 25);

    elseif nargin == 4

        if ischar(varargin{1})
            inSt = varargin{1};
            inInt = min(length(data), 25);
        elseif isnumeric(varargin{1})
            inInt = varargin{1};
            inSt = 'max';
        else
            disp('Warning: The fourth input must be either be a choice (max/min) or number of images to display. \n Setting default values...');
            inSt = 'max';
            inInt = min(length(data), 25); 
        end

    elseif nargin >= 5

        if ischar(varargin{1}) & isnumeric(varargin{2})
            inSt = varargin{1};
            inInt = varargin{2}; 
        else
            disp('Warning: The optional inputs must be a choice (max/min) and the number of images to display. \n Setting default values...');
            inSt = 'max';
            inInt = min(length(data), 25); 
        end

        if nargin == 6
            im_sel = varargin{3}; 
            cnt = 0;
        end

        if nargin > 6
            disp('Warning: Ignoring additional inputs.');
        end
    
    end

        
    if strcmp(inSt,'max')
        [~, im_idx] = sort(data, 'descend');
    elseif strcmp(inSt,'min')
        [~, im_idx] = sort(data, 'ascend');
    else 
        error('Expected input: max/min')
    end

    % Load and extract image filenames
    im_list = extractfield(dir(im_path),'name')';

    im_idx = shared_idx(im_idx);

    figure
    tiledlayout(ceil(length(shared_idx)/5),5)
    hold on

    for ii = 1 : inInt

        im_fname = im_list( contains( im_list, sprintf('shared%04d_nsd',im_idx(ii))));

        if isempty(im_fname)
            fprintf('Warning: Image shared%04d not found!',im_idx(ii))
        else
            nexttile
            if exist('cnt', 'var') && ismember(im_idx(ii), im_sel)
                imshowWithBorder(imread(fullfile(im_path, im_fname{1}))) 
                cnt = cnt+1;
            else
                imshow(imread(fullfile(im_path, im_fname{1})))
            end

        end

    end

    if exist('cnt', 'var')
        fprintf('\n Count selected images: %d \n', cnt);
    end

end

%%

function imshowWithBorder(img)

% Define border size
borderSize = 20; % Adjust as needed

% Create a red border
redBorder = uint8(255 * ones(size(img, 1) + 2 * borderSize, ...
                              size(img, 2) + 2 * borderSize, ...
                              size(img, 3))); % Create a red image of appropriate size
redBorder(:, :, 1) = 255; % Set red channel to maximum
redBorder(:, :, 2) = 0;   % Set green channel to zero
redBorder(:, :, 3) = 0;   % Set blue channel to zero

% Place the original image in the center of the red border
redBorder(borderSize + (1:size(img, 1)), ...
          borderSize + (1:size(img, 2)), :) = img;

% Display the image with red border
imshow(redBorder);

end