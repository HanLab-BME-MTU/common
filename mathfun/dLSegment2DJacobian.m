function [dFdXc dFdYc dFdA dFds dFdl dFdt] = ...
    dLSegment2DJacobian(xRange, yRange, xC, yC, A, sigmaPSF, l, theta)
% Partial derivatives of a 2-dimensional diffraction-limited segment.
% [dFdXc dFdXy dFdA dFds dFdl dFdt] =
%           dLSegment2D_dFdA(xRange, yRange, xC, xC, A, sigmaPSF, l, theta)
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
% dFdXc dFdXy dFdA dFds dFdl dFdt are NxM matrices where N = numel(xRange)
% and M = numel(yRange). 
%
% Sylvain Berlemont, 2010

ct = cos(theta);
st = sin(theta);

l = l / 2;

c0 = sqrt(2) * sigmaPSF;
c = A / (2 * erf(l / c0));

[X Y] = meshgrid(xRange, yRange);

X = X - xC;
Y = Y - yC;

dFdXc = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));

dFdYc = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));

dFdA = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));

dFds = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));

dFdl = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));

dFdt = c * exp(-((st * X - ct * Y) / c0).^2) .* ...
    (erf((l - ct * X - st * Y) / c0) + ...
    erf((l + ct * X + st * Y) / c0));
