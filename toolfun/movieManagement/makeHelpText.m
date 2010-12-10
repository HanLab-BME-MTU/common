function textOut = makeHelpText(description,paramList,paramDescriptions)
%MAKEHELPTEXT compiles the  input description and parameter description
%into a nicely formatted help text.
% 
% Input:
% 
%   description - A character array containing a general description of the
%   function and what it does.
% 
%   paramList - A cell array of character arrays containing a list of
%   parameters to create help sections for
% 
%   paramDescription - A cell array of character arrays containing
%   descriptions for each of the parameters in paramList
% 
%
% Output:
%
%   textOut - The compiled help text.
%
% Hunter Elliott
% 6/2010
%

if nargin < 1 || isempty(description)
    error('You must input a description!')
end

if nargin < 2
    paramList = {};
end

if nargin < 3
    paramDescriptions = {};
end

textOut = ['Process Description: \n\n',...
           description '\n\n'];
       
nParam = numel(paramList);

if numel(paramDescriptions) ~= nParam
    error('You must input a description for each parameter in the paramList!')
end

if nParam > 0
    
    textOut = [textOut '\nParameter Descriptions: \n'];

    for j = 1:nParam

       textOut = [textOut '\n' paramList{j} ':\n\n' paramDescriptions{j} '\n'];

    end
   
    textOut = [textOut '\n\n\n\n\n\n'];

end
