function listOfFiles = searchFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode)
%searchFiles is an utility to search for files containing a specific string
%
%SYNOPSIS listOfFiles = searchFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode)
%
%INPUT    includeString: string contained in the filenames you are looking for
%         excludeString (opt): string not contained in the filenames you are looking for
%         directory (opt): directory to search. if empty, current directory is searched (default)
%                               if 'ask', program asks for directory
%         includeSubDirectories (opt): whether to search subdirectories or not (0/{1})
%         selectionMode (opt): which file(s) to select if there are several files matching the
%                               includeString within one directory
%                              {'all'}-all; 'new'-newest; 'old'-oldest; 'GUI'-open GUI to select one file
%
%OUTPUT  listOfFiles: cell array with {[fileName] [directoryName]}
%
%c: 7-03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---test input---
if nargin<1|isempty(includeString)
    error('not enough input arguments or empty includeString')
end

%includeString
if ~isstr(includeString)
    error('includeString has to be a string')
end

%excludeString
if nargin>1&~isempty(excludeString)
    if ~isstr(excludeString)
        error('excludeString has to be a string')
    end
else
    excludeString = [];
end

%directory (ask if necessary)
if nargin>2
    if isempty(directory)
        directory = pwd;
    elseif strcmp(directory,'ask')
        directory = uigetdir(pwd,'select a directory to search');
        if isempty(directory)
            error('searchFiles aborted by user')
        end
    elseif ~isdir(directory)
        error([directory,' is not a valid directory!'])
    end
end
if nargin<3|isempty(directory)
    directory = pwd;
end

%includeSubDirectories
if nargin>3&~isempty(includeSubDirectories)
    if includeSubDirectories ~= 1
        includeSubDirectories = 0;
    end
else
    includeSubDirectories = 1;
end
if nargin<4|isempty(includeSubDirectories)
    includeSubDirectories = 1;
end

%selectionMode
if nargin<5|isempty(selectionMode)
    selectionMode = 'all';
else
    if ~(strcmp(selectionMode,'new')+strcmp(selectionMode,'old')+strcmp(selectionMode,'GUI')+strcmp(selectionMode,'all'))
        error('wrong selectionMode')
    end
end

%---end test input---


%---collect files---

%look for directories and wanted files in current directory, then check all subdirs

%init variables
dirs2check{1} = directory;
listOfFiles = [];
oldDir = pwd;

while ~isempty(dirs2check)
    %init/empty var
    pathCell = {};
    
    %open dir & store currentDir
    cd(dirs2check{1});
    currentDir = pwd;
    
    %list all files of current directory
    currentDirList = dir;
    
    %look for subdirectories and add to dirs2check
    isDirList = cat(1,currentDirList.isdir);
    if length(isDirList)>2
        subDirIdx = find(isDirList(3:end))+2;
        for i=1:length(subDirIdx)
            dirs2check = [dirs2check; {[currentDir,filesep,currentDirList(subDirIdx(i)).name]}];
        end
    end %if length(isDirList)>2
    
    %look for files in current directory and store them
    newFiles = chooseFile(includeString,[],selectionMode,excludeString);
    if ~isempty(newFiles)
        if ~iscell(newFiles)
            newFiles = cellstr(newFiles);
        end
        [pathCell{1:size(newFiles,1)}] = deal(currentDir);
        listOfFiles = [listOfFiles; newFiles, pathCell'];
    end %if ~isempty(newFiles)
    
    %eliminate top directory
    dirs2check(1)=[];
    
    %check wheter we want to look at subDirs
    if ~includeSubDirectories
        dirs2check = [];
    end
    
    
end %while ~isempty(dirs2check)

%---end collect files---

cd(oldDir)