function [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigma,l,theta,imSize)
% [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigma,l,theta,imSize)
%
% Compute the finite support of a diffraction-limited 2D segment given
% its parameters.
%
% parameters:
% (xC,yC)            center of the segment
%
% sigma              half width of segment
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% imSize             image size
%
% output:
% (xRange, yRange)   2 vectors representing the 2-dimensional support of
%                    the segment in the image domain 1:size(imSize,1) x
%                    1:size(imSize,2).
%
% nzIdx              linear indices where pixel value is not zero. These
%                    indices are local and are intended to be passed to
%                    subResSegment2D() and subResSegment2DJacobian()
%                    functions.
%
% Sylvain Berlemont, 2010

% half length of the segment
l2 = l / 2;

% Hypothenuse length, corresponding to the half-length diagonal of a
% 2*(L2+d) long by 2d wide rectangle surrounding the segment.
d = 3 * sigma;
lh = sqrt(d.^2 + (l2 + d).^2);

% Angle between a rectangle border and a diagonal of the rectangle
at = atan(d ./ (l2 + d));

s1 = [1 1 -1 -1];
s2 = [1 -1 1 -1];

% xy-coordinates of the 4 rectangle's corners.
x = xC + s1 .* cos(theta + s2 * at) * lh;
y = yC + s1 .* sin(theta + s2 * at) * lh;

% truncate numbers towards zero with 10 decimals.
x = fix(x * 1e10) * 1e-10;
y = fix(y * 1e10) * 1e-10;

xMin = min(floor(x));
xMax = max(ceil(x));
yMin = min(floor(y));
yMax = max(ceil(y));

xRange = max(xMin,1):min(xMax,imSize(2));
yRange = max(yMin,1):min(yMax,imSize(1));

% % Faster solution (maybe...)
% [X,Y] = meshgrid(xRange,yRange);
% ct = cos(theta);
% st = sin(theta);
% D1 = abs(ct * (xC - X) + st * (Y - yC)); % Wrong !
% D2 = abs(st * (xC - X) + ct * (Y - yC)); % Wrong !
% nzIdx = find(D1 <= d & D2 <= l2 + d);

% TODO: replace poly2mask call by 4 line intersections.
nzIdx = find(poly2mask(x-xRange(1)+1,y-yRange(1)+1,length(yRange),length(xRange)));
