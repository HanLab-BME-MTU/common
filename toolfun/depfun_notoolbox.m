function depList = depfun_notoolbox(funList)
%DEPFUN_NOTOOLBOX find dependencies of the input m-files, excluding toolbox and matlab built-in functions
% 
% depList = depfun_notoolbox(funList)
% 
% Returns a cell array with the paths of every m file which the input list
% of m files calls / depends on, including object classes, but excluding
% those functions which are provided by mathworks - toolbox functions,
% built-in matlab functions etc. The matlab function depfun will return
% toolboxes in a recursive search, so this function is used to find
% non-toolbox functions which are required.
% 
% Input:
% 
%   funList - a cell array of character strings containing the name(s) of
%   the files to check the dependencies of.
% 
% Output:
% 
%   depList - A cell array of character strings containing the full path
%   and file name of all the m-files the input funlist depends on.
%
% Hunter Elliott 6/2010
% Sebastien Besson July 2011

% Get the initial set of file dependencies
depList=cellfun(@which,funList(:),'UniformOutput',false);

while true
    %Find dependencies of current list
    newFiles = depfun(depList{:},'-toponly','-quiet');
        
    %Remove all the toolbox entries
    newFiles = newFiles(cellfun(@(x)(isempty(regexp(x,'toolbox','once'))),newFiles));    
    
    % Update the dependencies file list
    nFilesOld=numel(depList);
    depList = unique(vertcat(depList,newFiles));
    
    % Break if no new file or the same set of files is found
    if isempty(newFiles) || nFilesOld==numel(depList), return; end
end