function [depList, toolboxes] = getFunDependencies(funList, varargin)
%GETFUNDEPENDENCIES list dependencies and required toolboxes of the input m-files
% 
% SYNOPSIS depList = getFunDependencies(funList)
%          [depList,toolboxes] = getFunDependencies(funList)
%          [depList,toolboxes] = getFunDependencies(funList, excludefilter)
% 
% Returns a cell array with the paths of every m file which the input list
% of m files calls / depends on, including object classes, but excluding
% those functions which are provided by mathworks - toolbox functions,
% built-in matlab functions etc. The matlab function 
% matlab.codetools.requiredFilesAndProducts will return
% toolboxes in a recursive search, so this function is used to find
% non-toolbox functions which are required.
% Additionally returns the matlab toolboxes required to run these dependencies. 
% 
% Input:
% 
%   funList - a string or a cell array of strings containing the names of
%   the files whose dependencies should be listed.
% 
%   excludefilters - a string or cell array of regular expressions to
%   exclude from the dependency search. Default: empty cell array.
%
%   allMex - a logical. If true, then all versions of any  mex binaries
%   will be returned (any files with *.mex* in the filename, so long as
%   they are in the same folder as the binary for the current platform).
%   Default is false (only mex binaries for current platform are returned).
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
% Mark Kittisopikul, February 2016 depfun is now deprecated in lieu of
%                           matlab.codetools.requiredFilesAndProducts
% Based on depfun_notoolbox.m and toolboxesUsed.m

% Input check
if ischar(funList), funList={funList}; end
ip = inputParser;
ip.addRequired('funList', @iscell);
ip.addOptional('excludefilters', {}, @(x) ischar(x) || iscell(x));
ip.addOptional('allMex', false, @islogical);
ip.parse(funList, varargin{:});
if ischar(ip.Results.excludefilters)
    excludefilters = {'toolbox', ip.Results.excludefilters};
else
    excludefilters = [{'toolbox'},ip.Results.excludefilters{:}];
end

% Get the initial set of file dependencies
filesList=cellfun(@which,funList(:),'UniformOutput',false);
depList={};

% Construct anonymous functions for filtering files
filter = @(f, l) cellfun(@(x)(isempty(regexp(x,f,'once'))), l);

while true
    %Find dependencies of current list
    newFiles = matlab.codetools.requiredFilesAndProducts(filesList,'toponly');
        
    % Filter new found files using regular expression
    for i = 1:numel(excludefilters)
        newFiles = newFiles(filter(excludefilters{i}, newFiles));
    end
    
    % Update the dependencies file list
    nFiles=numel(depList);
    filesList = setdiff(newFiles,depList);
    depList = unique(vertcat(depList(:),newFiles(:)));
    
    % Break if no new file or the same set of files is found
    if isempty(newFiles) || nFiles==numel(depList), break; end
end

% Matlab toolbox files are named '*/toolbox/name_of_toolbox/*'
[~,productList] = matlab.codetools.requiredFilesAndProducts(depList(:),'toponly');
toolboxes = {productList.Name}';

if ip.Results.allMex
    %Find any mex binaries.
    mexInd = cellfun(@(x)(strfind(x,'.mex')),depList,'Unif',0);
    isMex = find(~cellfun(@isempty,mexInd))';
    mexFuns = {};
    for j = isMex
        %And find any other mex files in same directory
        newBins = dir([depList{j}(1:mexInd{j}+3) '*']);
        mexFuns = [mexFuns arrayfun(@(x)(which(newBins(x).name)),1:numel(newBins),'Unif',0)];               
    end
    depList(isMex) = [];%Remove duplicates
    depList = vertcat(depList,mexFuns(:));%Add all binaries
end
