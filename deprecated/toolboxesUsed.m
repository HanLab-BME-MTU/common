function toolBoxes = toolboxesUsed(funList)
%TOOLBOXESUSED - returns a list of all the MATLAB toolboxes dependend on by specified functions 
% 
% toolBoxes = toolboxesUsed(funList)
% 
% This function returns a list of all the MATLAB toolboxes required by the
% input list of functions. Note that the returned toolbox names are based
% on the folder the toolbox functions are stored in, not on the standard
% toolbox name.
% 
% Input:
% 
%   funList - A cell-array of character arrays containing the names of the
%   functions to check toolbox dependency for.
% 
% 
% Output:
% 
%   toolBoxes - A cell-array of character arrays containing the name(s) of
%   the toolbox(es) used by the list of input functions funList.
% 
% 
% Hunter Elliott
% 12/2010
%

fprintf(2, ['Warning: ''' mfilename ''' is deprecated and should no longer be used. Use getFunDependencies instead.\n']);

if nargin < 1 || isempty(funList)    
    error('Come on, what the fuck? You have to input SOMETHING!')
elseif strcmpi('SOMETHING',funList)
    error('Nice try, Smartass.')
elseif ~iscell(funList)
    error('The funList input must be a cell-array!')
end

disp('Determining function dependencies. Please be patient, this may take some time...')

%Note: This is probably about the slowest way possible to do this. Feel
%free to re-write it.

keepGoin = true;

toolBoxes = {};
funList = funList(:);
nFilesNew = numel(funList);
firstTime = true;

while keepGoin
    
    %Find dependencies of current list
    newFiles = depfun(funList{:},'-toponly','-quiet');

    %Toolbox files should be named '*/toolbox/name_of_toolbox/*'
    toolboxToken =[filesep 'toolbox' filesep '(\w+).*?' filesep];
    foundTokens=cellfun(@(x)regexp(x,toolboxToken,'tokens','once'),newFiles,'UniformOutput',false);

    %Store any toolboxes from this round
    toolBoxes= vertcat(toolBoxes,unique([foundTokens{:}])');
    
    %Remove the toolboxes so they won't be searched next round
    isToolBox=~cellfun(@isempty,foundTokens);
    newFiles(isToolBox) = [];    
       
    if ~isempty(newFiles)
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


toolBoxes = unique(toolBoxes);

%Remove the "toolboxes" that come with MATLAB by default
builtIn = {'matlab','local','compiler','control'};
isBuiltIn = cellfun(@(x)(any(strcmp(x,builtIn))),toolBoxes);
toolBoxes = toolBoxes(~isBuiltIn);

