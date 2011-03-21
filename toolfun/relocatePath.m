function newPath = relocatePath(oldPath,oldRootDir,newRootDir)
% RELOCATEPATH generates the name of a path which root directory has been
% modified
%
% 
% Input:
% 
%   oldpath - A string containing the path which should be replaced
% 
%   oldrootdir - A string containing the root directory of the path which
%   will be substituted by this function. This can be any portion of the
%   oldpath.
% 
%   newrootdir - A string containing the new root directory belonging or
%   not to the same OS.
% 
% 
% Output: 
%
%   newpath - A string containing the name of the path in the new root
%   directory. Should be the empty string if error occurs.
%
% Sebastien Besson, 03/2011
%

% Check the old root directory is contained within the old path
newPath='';
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
newPath=regexprep(oldPath,regexptranslate('escape',oldRootDir),regexptranslate('escape',newRootDir));
newPath=regexprep(newPath,oldFilesep,newFilesep);

end