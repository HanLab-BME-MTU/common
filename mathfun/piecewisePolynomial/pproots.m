function [ r ] = pproots( pp )
%pproots Find the roots of a piecewise polynomial
%
% INPUT
% pp - a piecewise polynomial, such as from mkpp or spline
%
% OUTPUT
% r - roots of the piecewise polynomial
%     If pp.dim == 1, then a column vector of the roots
%     If pp.dim ~= 1, then cell array of length pp.dim, each a column
%     vector of the roots
%
% See also roots, mkpp, spline
%
% Mark Kittisopikul, Ph.D.
% Jaqaman Lab
% UT Southwestern
% February 2017

interval_size = diff(pp.breaks);
if(pp.dim ~= 1)
    interval_size = repmat(interval_size,pp.dim,1);
end

r = NaN(size(pp.coefs)-[0 1]);
order = size(r,2);
coefs = pp.coefs;

for ii=1:size(coefs,1)
    % TODO: accelerate roots for small polynomials?
    temp = roots(coefs(ii,:));
    % Invalidate roots outside of interval breaks
    temp(temp < 0) = NaN;
    temp(temp > interval_size(ii)) = NaN;
    temp(end+1:order) = NaN;
    % Store
%     r(ii,1:length(temp)) = temp.';   
    r(ii,:) = temp.';
end

if(pp.dim == 1)
    % Add the left side value of breaks to root value
    r = bsxfun(@plus,r,pp.breaks(1:end-1).');
    % Only return valid roots
    r = r(~isnan(r));
else
    % Add the left side value of breaks to root value
    r = reshape(r,pp.dim,pp.pieces,pp.order-1);
    r = bsxfun(@plus,r,repmat(pp.breaks(1:end-1),pp.dim,1));
    r = reshape(r,pp.dim,pp.pieces*(pp.order-1));
    r = num2cell(r,2);
    r = cellfun(@(r) r(~isnan(r)),r,'UniformOutput',false);
end


end

