function [missingList]=findMissingFiles(completeList)
%findMissingFiles checks for essential Matlab files which are not found in the defined paths
%
%SYNOPSIS missingList=findMissingFiles(completeList)
%
%INPUT    completeList: cell array of functions on which program 'programName' depends, generated with the
%           command 'list=depfun(programName);' on a machine where the program works
%
%OUTPUT   missingList: cell array of all function names that are not found in the matlab paths
%         
% note: with "for i=1:size(missingList,1);copyfile(missingList{i},'#folder#');end", you
%       can easily copy all missing files to the folder #folder#
%
%c: 1/03, Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test input
if ~iscell(completeList)
    error('sorry, wrong input. Generate cell array ''completeList'' with depfun');
end

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
    if ~isempty(fileSeparators)
        fileName=currentFile(fileSeparators(end)+1:end);
        
        %search paths for fileName
        isItThere=which(fileName, '-all');
        
        %write fileName into output if not there
        if isempty(isItThere)
            missingList{k,1}=currentFile;
            k=k+1;           
        end
    end
end
missingList=sort(missingList);

