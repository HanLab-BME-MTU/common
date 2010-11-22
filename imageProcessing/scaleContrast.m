function out = scaleContrast(in, rangeIn, rangeOut)

if nargin<2 
    rangeIn = [min(in(:)) max(in(:))];
    rangeOut = [0 255];
end
if nargin<3 || isempty(rangeOut)
    rangeOut = [0 255];
end
if isempty(rangeIn)
    rangeIn = [min(in(:)) max(in(:))]; 
end

if rangeIn(2)-rangeIn(1) ~= 0
    out = (in-rangeIn(1)) / (rangeIn(2)-rangeIn(1)) * (rangeOut(2)-rangeOut(1)) + rangeOut(1);
else
    out = zeros(size(in));
end