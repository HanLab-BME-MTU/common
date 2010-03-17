function dFdYc = dLSegment2D_dFdYc(xRange, yRange, xC, yC, A, sigmaPSF, l, theta)
% Partial derivative of a 2-dimensional diffraction-limited degment
% function in function of yC.
% dFdYc = dLSegment2D_dFdYc(xRange, yRange, xC, xC, A, sigmaPSF, l, theta);
%
% parameters:
% (xRange, yRange)   2 vectors representing the 2-dimensional domain (e.g.
%                    xRange = -10:.1:10, yRange = -5:.1:5
%
% (xC,yC)            center of the segment
%
% A                  amplitude of the segment
%
% sigmaPSF           half width of the gaussian PSF model.
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% output:
% dFdYc is a NxM matrix where N = numel(xRange) and M = numel(yRange).
%
% Sylvain Berlemont, 2009

ct = cos(theta);
st = sin(theta);

l = l / 2;

c0 = sqrt(2) * sigmaPSF;
c = A / (2 * erf(l / c0));

[X Y] = meshgrid(xRange, yRange);

X = X - xC;
Y = Y - yC;

dFdYc = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));
