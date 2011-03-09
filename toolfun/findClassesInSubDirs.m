function classList = findClassesInSubDirs(topdir)

%FINDCLASSESINSUBDIRS find allclasses defined below a top directory
%subdirectories
% 
% classList = findClassesInSubDirs(topdir)
% 
% Finds recursively all the m-files in the subdirectories of an initial
% directory which name corresponds to a class. Returns a cell array with
% the name of all found classes.
%
% Input:
% 
%   topdir - a string containing the path of the initial directory.
% 
% Output:
% 
%   classList - a cell array of character strings containing the name(s) of
%   the found classes.
%
%
% Sebastien Besson
% 3/2011
%

dirs = getSubDirs(topdir);
dirs{end+1} = topdir;

classList = cell(1,length(dirs));
for k = 1:length(dirs)
    mFilesList = dir([dirs{k} filesep '*.m']);
    
    if ~isempty(mFilesList)
        ismFileClass = logical(arrayfun(@(x) exist(x.name(1:end-2),'class'), mFilesList));
        classList{k}=arrayfun(@(x) x.name(1:end-2),mFilesList(ismFileClass), 'UniformOutput', false);    
    end
end

classList = vertcat(classList{:});