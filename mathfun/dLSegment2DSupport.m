function [xRange, yRange] = dLSegment2DSupport(xC, yC, sigmaPSF, l, theta)
% Compute the finite support of 2D Diffraction-limited Segment Model given
% its parameters.
% [xRange yRange] = dLSegment2DSupport(xC, xC, sigmaPSF, l, theta)
%
% parameters:
% (xC,yC)            center of the segment
%
% sigmaPSF           half width of the gaussian PSF model
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% output:
% (xRange, yRange)   2 vectors representing the 2-dimensional domain
%
% Sylvain Berlemont, 2010

% half length of the segment
l2 = l / 2;

% Hypothenuse length, corresponding to the half-length diagonal of a
% 2*(L2+d) long by 2d wide rectangle surrounding the segment.
d = 4 * sigmaPSF;
lh = sqrt(d.^2 + (l2 + d).^2);

% Angle between a rectangle border and a diagonal of the rectangle
at = atan(d ./ (l2 + d));

s1 = [1 1 -1 -1];
s2 = [1 -1 1 -1];

% xy-coordinates of the 4 rectangle's corners.
x = xC + s1 .* cos(theta + s2 * at) * lh;
y = yC + s1 .* sin(theta + s2 * at) * lh;

xMin = min(floor(x));
xMax = max(ceil(x));
yMin = min(floor(y));
yMax = max(ceil(y));

if xMax < xMin
    xRange = xMax:xMin;
else
    xRange = xMin:xMax;
end

if yMax < yMin
    yRange = yMax:yMin;
else
    yRange = yMin:yMax;
end
