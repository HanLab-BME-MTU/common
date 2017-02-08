function [ ppderiv ] = ppdiff( pp )
%ppdiff Take the derivative of a piecewise polynomial, returning a new
%piecewise polynomial
%
% INPUT
% pp - a piecewise polynomial, such as from mkpp or spline
%
% OUTPUT
% ppderiv - a piecewise polynomial that is the derivative of pp
%
% See also diff, mkpp, spline
%
% Mark Kittisopikul, Ph.D.
% Jaqaman Lab
% UT Southwestern
% February 2017

% Multiply coefficients by exponents
dcoefs = bsxfun(@times,pp.coefs(:,1:end-1),pp.order-1:-1:1);

% If not singularly valued, then reshape to appropriate size for mkpp
if(pp.dim ~= 1)
    dcoefs = reshape(dcoefs,pp.dim,[],pp.order-1);
end

ppderiv = mkpp(pp.breaks,dcoefs,pp.dim);

end

