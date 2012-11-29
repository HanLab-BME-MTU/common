% [r, udata, sdata] = getMultiplicity(data) returns the occurrences/Multiplicity of the elements of 'data'
%
% Inputs:
%         data : n-dimensional input array
%
% Outputs: 
%          rep : # of occurrences for each element of 'data'
%        udata : sorted 1-D array of unique values in 'data'
%        sdata : sorted 1-D array of values in 'data'

% Francois Aguet, 03/02/2012 (modified on 10/29/2012)

function [rep, udata, sdata] = getMultiplicity(data)
sdata = sort(data(:))';
rep = diff([0 find([diff(sdata)~=0 1])]);
udata = unique(sdata);