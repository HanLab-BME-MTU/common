function [N,X] = histogram(data,factor)
% HISTOGRAM generates a histogram using the "optimal" number of bins
%
% If called with no output argument, histogram plots into the current axes
%
% SYNOPSIS [N,X] = histogram(data,factor)
%
% INPUT    data: vector of input data
%          factor: (opt) factor by which the bin-widths are multiplied
%
% OUTPUT   N   : number of points per bin
%          X   : center position of bins
%
% c: 2/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
data = data(:);

if nargin < 2 || isempty(factor)
    factor = 1;
end

nData = length(data);
if nData < 20
    warning('Less than 20 data points!')
    nBins = ceil(nData/4);
else


    % create bins with the optimal bin width
    % W = 2*(IQD)*N^(-1/3)
    interQuartileDist = diff(prctile(data,[25,75]));
    binLength = 2*interQuartileDist*length(data)^(-1/3)*factor;

    % number of bins: divide data range by binLength
    nBins = floor((max(data)-min(data))/binLength);

    if ~isfinite(nBins)
        nBins = length(unique(data));
    end

end

% histogram
if nargout > 0
    [N,X] = hist(data,nBins);
else
    hist(data,nBins);
end