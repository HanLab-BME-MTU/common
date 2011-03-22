function out = scaleContrast(in, rangeIn, rangeOut)

if nargin<2 || isempty(rangeIn)
    rangeIn = [min(in(:)) max(in(:))];
end
if nargin<3 || isempty(rangeOut)
    rangeOut = [0 255];
end

if rangeIn(2)-rangeIn(1) ~= 0
    out = (double(in)-rangeIn(1)) / diff(rangeIn) * diff(rangeOut) + rangeOut(1);
else
    out = zeros(size(in));
end