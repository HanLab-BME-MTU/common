% [r, udata, sdata] = getMultiplicityInt(data) returns the occurrences/Multiplicity of the elements of 'data'
%
% Inputs:
%         data : n-dimensional input array
%
% Outputs: 
%          rep : # of occurrences for each element of 'data'
%        udata : sorted 1-D array of unique values in 'data'
%        sdata : sorted 1-D array of values in 'data'
%
% Note: NaN/Inf elements in input data are ignored
%
% This is an optimization for getMultiplicity for arrays of integers in a 
% narrow range.
%
% Alternatives:
% 1) getMultiplicity(data)
% 2) histc(data,unique(data))
% 3) accumarray(data(:)',1)

% Mark Kittisopikul, 2014/08/11

function [rep, udata, sdata] = getMultiplicityInt(data)

if(~isinteger(data))
    warning(['getMultiplicityInt: Non-integer data type ' class(data) ' given']);
    data = data(isfinite(data));
end

data = data(:);
minValue = min(data);
if(minValue ~= 1)
    data = data + ( 1 - minValue );
end

%where the magic happens
rep = accumarray(data,1);


% if unique values requested, then return those with nonzero repeitions
if(nargout > 1)
    maxValue = max(data) - (1 - minValue);
    udata = minValue : maxValue;
    udata = udata(rep ~= 0);
end

if(nargout > 2)
    % do run length decompression like rude in File Exchange
    repcum = cumsum(rep(rep ~= 0));
    idx = data;
    idx(1) = 1;
    idx(2:end) = 0;
    idx( repcum(1:end-1) + 1) = 1;
    idx = cumsum(idx);
    sdata = udata(:,idx);
end

rep = rep(rep ~= 0)';

end