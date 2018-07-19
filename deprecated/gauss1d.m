function y = gauss1d(x,s)
%GAUSS1D return the gaussian bell curve for the values in x
%
% SYNOPSIS y = gauss1d(x,s)
%
%   where y = 1/sqrt(s*pi) * exp(-1/2 * (x/s)^2);

fprintf(2,['Warning: ''' mfilename ''' is deprecated. Use normpdf instead.\n']);


y = exp(-1/2 * (x/s).^2);
y = y ./(sqrt(2*pi)*s);