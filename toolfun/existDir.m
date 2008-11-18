function success = existDir(dirName)
%EXISTDIR checks whether a directory exists
%
% SYNOPSIS: success = existDir(dirName)
%
% INPUT dirName: string or cell array of strings containing possible
%                directoryNames (either relative to current directory or
%                absolute paths
%
% OUTPUT success: logical 1 (or 0) if dirName is a directory (or not)
%
% REMARKS ExistDir performs the same function as exist(dirName,'dir'), but
%         it is much faster.
%
% created with MATLAB ver.: 7.7.0.471 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 17-Nov-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test input
if nargin < 1
    error('please supply directory name to existDir')
end

% check for cell array
if iscell(dirName)
    success = cellfun(@existDir,dirName);
else
    % try to cd to the directory. If it works, it's a directory!
    success = false;
    try %#ok<TRYNC>
        oldDir = cd(dirName);
        cd(oldDir);
        % if we arrive here without problem, it's a success
        success = true;
    % catch
        % it's not an accessible directory
    end
end