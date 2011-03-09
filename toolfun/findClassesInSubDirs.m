%fileList = findClassesInSubDirs(topdir)
% Finds all user-defined classes in all subdirectories of 'topdir'.

% Sebastien Besson, March 9, 2011

function classList = findCla ssesInSubDirs(topdir)

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