function [missingList,completeList]=findMissingFiles(completeList,excludeDirs,verbose)
%findMissingFiles checks for essential Matlab files which are not found in the defined paths
%
%SYNOPSIS missingList=findMissingFiles(completeList,excludeDirs)
%
%INPUT    completeList: cell array of functions on which program 'programName' depends, generated with the
%           command 'list=depfun(programName);' on a machine where the program works
%           If completeList is empty, a dialog box allows you to search for
%           files.
%         excludeDirs : directories to exclude in the path search. The
%           match is case-sensitive, and the drive name is lower case.
%           If you have specified excludeDirs in the subfunction, you can
%           instead supply a numeric selection.
%         verbose : verbosity of fdep. Optional. Default: false
%
%OUTPUT   missingList: cell array of all function names that are not found in the matlab paths
%         completeList: complete list of file dependencies
%
% note: with "for i=1:size(missingList,1);copyfile(missingList{i},'#folder#');end", you
%       can easily copy all missing files to the folder #folder#
%
%c: 1/03, Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(verbose)
    verbose = false;
end
if verbose
    verboseString = '';
else
    verboseString = '-q';
end

%test input
if isempty(completeList)
    bait = uipickfiles('REFilter','\.m$');
    completeList = {};
    try
        for b = 1:length(bait)
            fList = fdep(bait{b},verboseString);
            completeList = unique([completeList;fList.froot]);
        end
        % fix file separator orientation
        for c = 1:length(completeList)
            completeList{c} = strrep(completeList{c},'/',filesep); %#ok<AGROW>
        end
    catch err
        if strcmp(err.identifier,'MATLAB:UndefinedFunction')
            disp('It is suggested to use fdep (from FileExchange) instead of depfun for better performance')
            completeList = depfun(bait{:});
        else
            rethrow(err)
        end
    end
end
if ~iscell(completeList)
    error('sorry, wrong input. Generate cell array ''completeList'' with depfun');
end
if nargin < 2 || isempty(excludeDirs)
    excludeDirs = [];
elseif isnumeric(excludeDirs)
    excludeDirs = defaultExcludes(excludeDirs);
elseif ~iscell(excludeDirs)
    excludeDirs = {excludeDirs};
end
nExclude = numel(excludeDirs);


%initialize variables
k=1;
missingList{1,1}='';

%loop through completeList to build missingList
for i=1:size(completeList,1)
    %get first string
    currentFile=completeList{i};
    %extract filename
    fileSeparators=findstr(currentFile,filesep);
    %fileName: all chars after the last file separator. If no file separator, it is no file
    % also, check that the file is not on the matlab path
    if isempty(strmatch(matlabroot,currentFile)) &&  ~isempty(fileSeparators)
        fileName=currentFile(fileSeparators(end)+1:end);
        
        %search paths for fileName
        isItThere=which(fileName, '-all');
        
        % check whether to exlcude dirs
        iter = 1;
        while ~isempty(isItThere) && iter <= nExclude
            % match all or beginning
            idx = strmatch(excludeDirs{iter},isItThere);
            if ~isempty(idx)
                isItThere(idx) = [];
            end
            iter = iter + 1;
        end
        
        %write fileName into output if not there
        if isempty(isItThere)
            missingList{k,1}=currentFile;
            k=k+1;
            % check whether there is a fig file with the same name
            [pn,fn] = fileparts(currentFile);
            figFile = fullfile(pn,fn,'.fig');
            if exist(figFile,'file')
                missingList{k,1} = figFile;
                k = k + 1;
            end
        end
    end
end
missingList=sort(missingList);

% loop through missingList. If no extension, add .m
for m = 1:length(missingList)
    [p,f,e]=fileparts(missingList{m});
    if isempty(e)
        missingList{m} = [missingList{m},'.m'];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
function excludeDirs = defaultExcludes(selection)
% defaultExcludes is a place to store default exclude directory lists

switch selection
    case 1
        excludeDirs = {'c:\data\jonas\matlab\common-tsri';...
            'c:\data\jonas\matlab\extern-tsri';...
            'c:\data\jonas\matlab\newFunctions';...
            'c:\data\jonas\matlab\mdxMisc';...
            'c:\data\jonas\matlab\chromdyn-tsri'};
    case 2
        excludeDirs = {'/Users/jonas/matlab'};
    case 3
        excludeDirs = {'/Users/jonas/matlab/common-tsri';...
            '/Users/jonas/matlab/extern-tsri';...
            '/Users/jonas/matlab/chromdyn-tsri';...
            '/Users/jonas/matlab/lccbCommon-hms';...
            '/Users/jonas/matlab/lccbExtern-hms'};
    case 4
        excludeDirs = {'c:\data\matlab\lccbCommon-iric';...
            'c:\data\matlab\lccbMaki-iric';...
            'c:\data\matlab\lccbChromdyn-iric';...
            'c:\data\matlab\qu_sf';...
            'c:\data\matlab\newFunctions';...
            'c:\data\matlab\tanTrao';...
            'c:\data\matlab\daniel';...
            'c:\data\matlab\mdxMisc'};
    case 5 
        excludeDirs = {'c:\data\matlab'};
    otherwise
        excludeDirs = [];
end

