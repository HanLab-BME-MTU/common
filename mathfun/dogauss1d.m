function y = dogauss1d(x,s)
%DOGAUSS1D returns the 1st derivative of the gaussian bell curve for the values in x
%
% SYNOPSIS y = dogauss1d(x,s)
%
%   where x = 1/sqrt(s*pi) * (-x/s^2) * exp(-1/2 * (x/s)^2);
% 
y = exp(-1/2 * (x/s).^2) ./ (sqrt(2*pi)*s);
y = -y.*x/(s^2);