function [ maxima, maxima_value, minima, minima_value, other, other_value ] = ppextrema( pp )
%ppextrema Find the extrema of a piecewise polynomial
%
% INPUT
% pp - piecewise polynomial
%
% OUTPUT
% maxima - local maxima, not including boundary points
% maxima_value - value of pp at maxima
% minima - local minima, not including boundary points
% minima_value - value of pp at minima
% other - critical points where first and second derivative is zero
% other_value - value of pp at other
%
% See also mkpp, spline, ppdiff, pproots
%
% Mark Kittisopikul, February 2017
% Jaqaman Lab
% UT Southwestern

TOL = eps('single');

% First derivative
ppd = ppdiff(pp);
% Roots of first derivative
ppd_r = pproots(ppd);
if(isreal(pp.coefs))
    ppd_r = ppd_r(abs(imag(ppd_r)) < TOL);
end

% Second derivative
ppdd = ppdiff(ppd);
ddval = ppval(ppdd,ppd_r);

% Maxima
maxima = ppd_r(ddval < 0);

if(nargout > 1)
    maxima_value = ppval(pp,maxima);
end

if(nargout > 2)
    minima = ppd_r(ddval > 0);
end

if(nargout > 3)
    minima_value = ppval(pp,minima);
end

if(nargout > 4)
    other = ppd_r(ddval == 0);
end

if(nargout > 5)
    other_value = ppval(pp,other);
end

end

