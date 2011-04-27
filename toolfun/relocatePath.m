function newPath = relocatePath(oldPath,oldRootDir,newRootDir)
% RELOCATEPATH relocates all paths in a given input
%
% 
% Input:
% 
%   oldpath - The element to be relocated by the function. The function
%   currently supports three types of input:
%       1 - a string
%       2 - a structure or structure array 
%       3 - a cell array
% 
%   oldrootdir - A string containing the root directory of the path which
%   will be substituted by this function. This can be any portion of the
%   oldpath.
% 
%   newrootdir - A string containing the new root directory (belonging or
%   not to the same OS).
% 
% 
% Output: 
%
%   newpath - The function output depends of the type of input
%       1 - a string containing the name of the relocated path 
%       2 - a structure or structure array where all string fields have
%       been relocated
%       3 - a cell array where all strings have been relocated 
%
% Sebastien Besson, 03/2011
%

% Check argument number and type
if nargin < 3  
    error('Three non-empty inputs required!')
end

if ~ischar(oldRootDir) || ~ischar(newRootDir)
    error('Inputs 2 and 3 must be character strings!')
end

% Return the input by default (no relocation)
newPath=oldPath;
if (isempty(oldRootDir) && isempty(newRootDir)), return; end

% Recursive function call if input is structure array or cell array
if isstruct(oldPath)    
    newPath = arrayfun(@(x) structfun(@(y) relocatePath(y,oldRootDir,newRootDir),x,'UniformOutput',false),oldPath);
    return
elseif iscell(oldPath)
    newPath=cellfun(@(x) relocatePath(x,oldRootDir,newRootDir),oldPath,'UniformOutput',false);
    return
end

% Check the old root directory is contained within the old path
if numel(oldPath)<numel(oldRootDir)
    return
elseif ~strcmp(oldPath(1:numel(oldRootDir)),oldRootDir);
    return;
end

% Get file separators of old and new root directories as regular
% expressions
oldFilesep=getFilesep(oldRootDir);
newFilesep=getFilesep(newRootDir);

% Remove ending separators in the paths
oldPath=regexprep(oldPath,[oldFilesep '$'],'');
oldRootDir=regexprep(oldRootDir,[oldFilesep '$'],'');
newRootDir=regexprep(newRootDir,[newFilesep '$'],'');

% Generate the new path and replace the wrong file separators
if isempty(oldRootDir)
    % In the case of mount relocation from /mnt to /home/xxx/mnt
    newPath = [newRootDir oldPath];
else
    newPath=regexprep(oldPath,regexptranslate('escape',oldRootDir),regexptranslate('escape',newRootDir));
end
newPath=regexprep(newPath,oldFilesep,newFilesep);

end