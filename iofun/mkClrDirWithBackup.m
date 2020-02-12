function [] = mkClrDirWithBackup(dirPath)
%mkClrDirWithBackup makes a folder using mkClrDir but backs up whatever
%therre is files in the existing file in dirPath
% 
% This is just a little function for creating / settin up output
% directories. It checks if a directory exists and if not, makes it. If the
% directories does exist and contains files, these are deleted.
% 
% Input:
% 
%   dirPath - the path to the directory to make/clear.
% 
% Sangyoon Han
% 11/2019

if nargin < 1 || isempty(dirPath)
    error('You must specify the directory to set up!')
end

if strcmp(dirPath(end),filesep)
    dirPath=[dirPath filesep];
end

if exist(dirPath,'dir')
    disp('Backing up the original data')
    ii = 1;
    backupFolder = [dirPath ' Backup ' num2str(ii)];
    while exist(backupFolder,'dir') && isempty(dir(backupFolder))
        backupFolder = [dirPath ' Backup ' num2str(ii)];
        ii=ii+1;
    end

    movefile(dirPath, backupFolder, 'f')
end
mkClrDir(dirPath);

end

