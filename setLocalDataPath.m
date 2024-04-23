function localDataPath = setLocalDataPath(varargin)
% function LocalDataPath = setLocalDataPath(varargin)
% Return the path to the root  directory and add paths in this repo
%
% input:  
%   personalDataPath: optional, set to 1 if adding personalDataPath
%
% when adding personalDataPath, the following function should be in the
% root of this repo:
%
% function localDataPath = personalDataPath()
%     'localDataPath = [/my/path/to/data];
%             
% function localDataPath = personalDataPath(1)
%     % to load path from personalDataPath
%
% personalDataPath is ignored in .gitignore
%
% dhermes, 2020, Multimodal Neuroimaging Lab


if isempty(varargin)
    rootPath = which('setLocalDataPath');
    RepoPath = fileparts(rootPath);

    % add path to functions
    addpath(genpath(RepoPath));

    % add localDataPath default
    localDataPath = fullfile(RepoPath,'data');
elseif ~isempty(varargin)
    % add path to data
    if varargin{1}==1 && exist('personalDataPath','file')
        localDataPath = personalDataPath();
        
    elseif varargin{1}==1 && ~exist('personalDataPath','file')
        sprintf(['add personalDataPath function to add your localDataPath:\n'...
            '\n'...
            'function localDataPath = personalDataPath()\n'...
            'localDataPath.input = [/my/path/to/data];\n'...
            'localDataPath.output = [/my/path/to/output];\n'...
            '\n'...
            'this function is ignored in .gitignore'])
        return
    end
    
    % add path to functions
    rootPath = which('setLocalDataPath');
    RepoPath = fileparts(rootPath);
    addpath(genpath(RepoPath));
end

return

