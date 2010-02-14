function moveSubsetOfFiles(dirSourceTarget,indices2move,filenameBase,...
    filenameExt,digits4Enum,renumber)
%MOVESUBSETOFFILES moves a subset of files from one directory to another, with the possibility of renumbering
%
%SYNPOSIS moveSubsetOfFiles(dirSourceTarget,indices2move,filenameBase,...
%    filenameExt,digits4Enum,renumber)
%
%INPUT  dirSourceTarget: Structure with fields:
%           .source : Directory to move files from.
%           .target : Directory to move file to.
%       indices2move: Row vector of file indices to move.
%       filenameBase: File name base (everything before the numbers).
%       filenameExt : File name extension (everything after the ".").
%       digits4Enum : Number of digits used to enumerate files.
%       renumber    : 1 to renumber the moved files, starting with 1 and
%                     using increments of 1; 0 otherwise. When files are
%                     renumbered, "new_" is inserted after the filenameBase
%                     and before the numbers.
%                     Optional. Default: 0.
%
%OUTPUT NONE
%
%Khuloud Jaqaman, February 2010

%% Input

%check whether correct number of input arguments was used
if nargin < 5
    disp('--moveSubsetOfFiles: Incorrect number of input arguments!');
    return
end

%get source and target directories
sourceDir = dirSourceTarget.source;
targetDir = dirSourceTarget.target;

%check whether to renumber files
if nargin < 6 || isempty(renumber)
    renumber = 0;
end

%% Action!

%based on renumbering choice ...
if renumber %if file names are renumbered
    
    %initialize index for renumbering
    iNew = 0;
    
    %go over files to be moved
    for iFile = indices2move
        
        %get file name
        file2move = [filenameBase sprintf(['%0' num2str(digits4Enum) 'i'],iFile) '.' filenameExt];
        
        %construct new file name
        iNew = iNew + 1;
        newName = [filenameBase 'new_' sprintf(['%0' num2str(digits4Enum) 'i'],iNew) '.' filenameExt];
        
        %move file
        movefile(fullfile(sourceDir,file2move),fullfile(targetDir,newName));
        
    end
        
else %if file names stay the same
    
    %go over files to be moved
    for iFile = indices2move
        
        %get file name
        file2move = [filenameBase sprintf(['%0' num2str(digits4Enum) 'i'],iFile) '.' filenameExt];
        
        %move file
        movefile(fullfile(sourceDir,file2move),targetDir);
        
    end
    
end

%% ~~~ the end ~~~
