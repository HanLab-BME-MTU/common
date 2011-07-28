function funList = depfun_notoolbox(funList)
%DEPFUN_NOTOOLBOX find dependencies of the input m-files, excluding toolbox and matlab built-in functions
% 
% funList = depfun_notoolbox(funList)
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
%   funList - A cell array of character strings containing the full path
%   and file name of all the m-files the input funlist depends on.
%
%
% Hunter Elliott
% 6/2010
%

fprintf(2, ['Warning: ''' mfilename ''' is deprecated and should no longer be used. Use getFunDependencies instead.\n']);

%Make sure funList is a column
funList = funList(:);
nFinit = numel(funList);
nFilesNew = nFinit;

keepGoin = true;
firstTime = true;

while keepGoin
    
    %Find dependencies of current list
    newFiles = depfun(funList{:},'-toponly','-quiet');
        
    
    %Remove all the toolbox entries
    newFiles = newFiles(cellfun(@(x)(isempty(regexp(x,'toolbox','ONCE'))),newFiles));    
    
    if ~isempty(newFiles) && ~firstTime        
        funList = vertcat(funList,newFiles);
        funList = unique(funList);        
    elseif firstTime
        firstTime = false;
        %Don't include the input list, as these will be duplicated
        funList = newFiles;
    else
        keepGoin = false;
    end
    
    nFilesOld = nFilesNew;
    nFilesNew = numel(funList);
    
    if nFilesOld == nFilesNew
        keepGoin = false;
    end
    
    
end