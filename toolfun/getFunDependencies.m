function [depList toolboxes] = getFunDependencies(funList)
%GETFUNDEPENDENCIES list dependencies and required toolboxes of the input m-files
% 
% SYNOPSIS depList = getFunDependencies(funList)
%          [depList,toolboxes] = getFunDependencies(funList)
% 
% Returns a cell array with the paths of every m file which the input list
% of m files calls / depends on, including object classes, but excluding
% those functions which are provided by mathworks - toolbox functions,
% built-in matlab functions etc. The matlab function depfun will return
% toolboxes in a recursive search, so this function is used to find
% non-toolbox functions which are required.
% Additionally returns the matlab toolboxes required to run these dependencies. 
% 
% Input:
% 
%   funList - a string or a cell array of strings containing the names of
%   the files whose dependencies should be listed.
% 
% Output:
% 
%   depList - A cell array of strings containing the full path of the
%   dependencies of the input, excluding the matlab functions.
%
%   toolBoxes - A cell array of character strings containing the name of
%   the matlab toolboxes which the depList depends on.

% Hunter Elliott,  June 2010
% Sebastien Besson, July 2011
% Based on depfun_notoolbox.m and toolboxesUsed.m

% Input check
if ischar(funList), funList={funList}; end
ip = inputParser;
ip.addRequired('funList',@iscell);
ip.parse(funList);

% Get the initial set of file dependencies
filesList=cellfun(@which,funList(:),'UniformOutput',false);
depList={};

while true
    %Find dependencies of current list
    newFiles = depfun(filesList{:},'-toponly','-quiet');
        
    %Remove all the toolbox entries
    newFiles = newFiles(cellfun(@(x)(isempty(regexp(x,'toolbox','once'))),newFiles));    
    
    % Update the dependencies file list
    nFiles=numel(depList);
    filesList = setdiff(newFiles,depList);
    depList = unique(vertcat(depList,newFiles));
    
    % Break if no new file or the same set of files is found
    if isempty(newFiles) || nFiles==numel(depList), break; end
end

% Matlab toolbox files are named '*/toolbox/name_of_toolbox/*'
allDepFiles = depfun(depList{:},'-toponly','-quiet');
toolboxToken = ['toolbox' regexptranslate('escape',filesep) '(\w+)' regexptranslate('escape',filesep)];
foundTokens=regexp(allDepFiles,toolboxToken,'tokens','once');
toolboxes= unique(vertcat(foundTokens{:}));

% Remove the "toolboxes" that come with MATLAB by default
builtInToolboxes = {'matlab','local','compiler','control'};
toolboxes = setdiff(toolboxes,builtInToolboxes);
