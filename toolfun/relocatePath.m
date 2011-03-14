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
if numel(oldPath)<numel(oldRootDir)
    return
elseif ~strcmp(oldPath(1:numel(oldRootDir)),oldRootDir);
    return;
end

% Identify old and new separators
oldFilesep=unique(regexp(oldRootDir,'/|\','match','once'));
newFilesep=unique(regexp(newRootDir,'/|\','match','once'));

% Remove ending separators in the paths
endingOldFilesepToken = [regexptranslate('escape',oldFilesep) '$'];
endingNewFilesepToken = [regexptranslate('escape',newFilesep) '$'];
oldPath=regexprep(oldPath,endingOldFilesepToken,'');
oldRootDir=regexprep(oldRootDir,endingOldFilesepToken,'');
newRootDir=regexprep(newRootDir,endingNewFilesepToken,'');

% Generate the new path and replace the extra file separator
newPath=regexprep(oldPath,regexptranslate('escape',oldRootDir),regexptranslate('escape',newRootDir));
newPath=regexprep(newPath,oldFilesep,newFilesep);

end
