function [uniqueEntries,numberOfOccurences] = countEntries(m)
%COUNTENTRIES returns all unique entries (sorted) in m and how many times the respective entries occured
%
%SYNOPSIS [uniqueEntries,numberOfOccurences] = countEntries(m)
%
%INPUT  m: any matrix (not cells or structs)
%
%OUTPUT uniqueEntries : unique(m)
%       numberOfOccurences : how many times the unique entries appear in m
%
%c: 11/03, jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test input
if iscell(m) | isstruct(m)
    error('cells and structs are not supportet as input');
end

%make m into a vector
m = m(:);

%get unique Entries
[uniqueEntries, dummy, uniqueIdx] = unique(m);
uniqueIdx = sort(uniqueIdx);

%uniqueIdx now has e.g. the form [1,1,1,2,3,3,3,4,4], so the function diff will
%return 1 wherever we change from one index to another
diffUniqueIdx = diff(uniqueIdx);
idxJump = find(diffUniqueIdx); %[3;4;7]

%if we now add 0 and the length of uniqueIdx to idxJump, another diff will return the
%number of equal indices we had in uniqueIdx.
idxJump = [0;idxJump;length(uniqueIdx)]; %[0;3;4;7;9]

numberOfOccurences = diff(idxJump); %[3;1;3;2]

