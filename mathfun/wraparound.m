function [wrappedValues, multiplier] = wraparound(values,interval)
%WRAPAROUND wraps values onto an interval [a, b)
%
% SYNOPSIS [wrappedValues, multiplier] = wraparound(values,interval)
%
% INPUT    values   : n-by-d array of values to be wrapped
%          interval : 2-by-1 o4 2-by-d array of intervals. The first row is
%                       the lower, the second row the upper limit of the
%                       interval. If a 2-by-1 interval is given with a
%                       n-by-d array of values, the same interval is
%                       applied for all dimensions
%
% OUTPUT   wrappedValues : n-by-d values as projected onto the interval
%          multiplier    : multiplier of interval to reconstitute original
%                           value
%
% REMARKS   The final value will lie within the interval [lower, upper).
%           It will never have th value of upper. To find out how
%           far away from the lower limit the value lies, you have to
%           subtract that lower limit yourself
%
% EXAMPLES  With an interval [0;10], 19 becomes [[9], [1]]
%           With an interval [2;4],  [9.2;-2.6] becomes [[3.2;3.4], [3;-3]]
%           With an interval [-180,180], [270;360] becomes [[-90;0], [1;1]]
%
% c: 7/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=====================
% TEST INPUT
%=====================

if nargin ~= 2 || isempty(values) || isempty(interval)
    error('You have to specify two non-empty input arguments!')
end

% get number of dimensions
sizeValues = size(values);
if length(sizeValues) > 2
    error('only 2D arrays allowed for values!')
end

% make sure interval is a 2-by-nDims array
interval = returnRightVector(interval, 2, 'r');
% fill up interval if necessary
if size(interval,2) == 1 && sizeValues(2) > 1
    interval = repmat(interval, [1, sizeValues(2)]);
end
if ~all(interval(1,:) < interval(2,:))
    error('lower limits have to be strictly smaller than upper limits!')
end

% make arrays filled with lower and upper limits
lower = repmat(interval(1,:),[sizeValues(1),1]);
upper = repmat(interval(2,:),[sizeValues(1),1]);

%=========================
% WRAP AROUND
%=========================

% Strategy: Transform interval and data to [0,1], then separate integer
% part as multiplier and transform the rest back onto the original interval
% With negative data, add 1 to the value
% shift values and upper limit by lower limit
wrappedValues = values - lower;
upper = upper - lower;

% find negative values
isNegative = wrappedValues < 0;

% divide by upper limit
wrappedValues = wrappedValues ./ upper;

% make multiplier. Floor rounds toward -inf, which works well with negative
% values
multiplier = floor(wrappedValues);

% now get remainder. Add 1 wherever we would get a negative wrapped value
wrappedValues = abs(wrappedValues) - multiplier + isNegative;

% transform wrappedValues back onto the original interval
wrappedValues = (wrappedValues .* upper) + lower;