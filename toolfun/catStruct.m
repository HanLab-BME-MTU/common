function out = catStruct(dim,structName)
%CATSTRUCT catenates values from multilevel structures
%
% SYNOPSIS out = catStruct(dim,structName)
%
% INPUT   dim        : dimension along which you want the struct to be catenated
%         structName : (string) full tree of the structure down to the field to be
%                      catenated, e.g. 'struct.fieldname.subfieldname'. 
%                      DO NOT USE COLONS OR RANGES
%
% OUTPUT  out        : array containtin the catenated values
%
% SAMPLE CALL : speed = catStruct(1,'wt30.individualStatistics.summary.antipolewardSpeed');
%
% currently, the software reassigns the output variable in every iteration of
% the loop. This can be slow. The function is relatively easily made faster
% though: preallocate out with size of field and levelsize. Two tests would
% be needed: first, go through the struct to find the first non-empty
% field, second, check for the size of out in the final eval and reassign
% if necessary.
% other possible improvements include support of ranges or indexLists for
% the fields (instead of doing full loops)
%
% c: 2/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=====================
% test input
%=====================

% nargin
if nargin ~=2 | isempty(dim) | isempty(structName)
    error('CATSTRUCT needs two non-empty input arguments!')
end

% structName
if ~isstr(structName)
    error('please input the structure as a string')
end

if ~strfind(structName,'.')
    error('no ''.'' found in structName. Make sure you input a structure')
end
if strfind(structName,'(') | strfind(structName,';') | strfind(structName,':') | strfind(structName,',') | strfind(structName,')')
    error('CATSTRUCT does not support ranges (yet)')
end

% dim is tested while catenating

%====================
% end test input
%====================



%=================================
% read structure info and catenate
%=================================

% parse structure name to find number of levels & number of entries. if
% only one level, we can use good old cat, otherwise, we will build a
% nested loop (horribile dictu!) to fill in an array

% '.' marks boundaries between levels
levelBreaks = strfind(structName,'.');

% a level is everything above the last field name
numberOfLevels = length(levelBreaks);

% load structure into mFile; make sure it gets the right name!!
topLevelName = structName(1:levelBreaks(1)-1);
eval([topLevelName '= evalin(''base'', [topLevelName '';'']);'])

% decide whether we can go easy or not
if numberOfLevels == 1
    % normal cat
    evalString = ['out = cat(' num2str(dim) ',' structName ');'];
else
    
    
    
    % adjust levelBreaks so that we can include the last filename
    levelBreaks = [levelBreaks length(structName)+1];
    
    % build nested loop
    startString = [];
    middleString = ['out = cat(' num2str(dim) ', out, '];
    endString   = []; 
    currentLevelName = topLevelName;
    % levelSize = 1; 
    
    for nLevel = 1:numberOfLevels
        
        % get next level name
        nextLevelName = [currentLevelName '(i' num2str(nLevel) ').' structName(levelBreaks(nLevel)+1:levelBreaks(nLevel+1)-1)];
        
        % % get number of elements of current level
        % levelSize = levelSize * eval(['numel(' currentLevelName ');']);
        
        % write loop
        startString = [startString 'for i' num2str(nLevel) '= 1:numel(' currentLevelName '),if ~isempty(' nextLevelName '),'];
        endString = ['end,end,' endString];
        
        % update currentLevelName
        currentLevelName = nextLevelName;
        
    end % for nLevel = 1:numberOfLevels
    
    % finish loop
    evalString = ['out = [];' startString middleString currentLevelName ');' endString];
end

% catenate structure

try
    eval(evalString);
catch
    rethrow(lasterr)
end

