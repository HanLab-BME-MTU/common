function F = dLSegment3D(X, Y, Z, xc, yc, zc, A, sigmaPSF, l, theta1, theta2)
% Diffraction-limited Segment Model
% F = dLSegment3D(X, Y, Z, xc, yc, zc, A, sigmaPSF, l, theta1, theta2)
%
% parameters:
% (X, Y, Z)    3 vectors representing the 2-dimensional domain (e.g X =
%              -10:.1:10, Y = -5:.1:5, Z = -5:.1:5
%
% (xc,yc,zc)   center of the segment
%
% A            amplitude of the segment
%
% sigmaPSF     half width of the gaussian PSF model.
%
% l            length of the segment
%
% theta1       orientation of the segment on xOy plan
%
% theta2       orientation of the segment on x0z plan
%
% output:
% F is a NxMxL matrix where N = numel(X), M = numel(Y) and L = numel(Z)
%
% Sylvain Berlemont, 2009

ct1 = cos(theta1);
ct2 = cos(theta2);
st1 = sin(theta1);
st2 = sin(theta2);

l = l / 2;

[X Y Z] = meshgrid(X,Y,Z);

X = X - xc;
Y = Y - yc;
Z = Z - zc;

c0 = sqrt(2) * sigmaPSF;
c = A / (2 * erf(l / c0));

tmp1 = l + Z * ct1 + X * ct2 * st1 + Y * st1 * st2;
tmp2 = l - Z * ct1 + X * ct2 * st1 - Y * st1 * st2;

F = c * exp(-(X.^2 + Y.^2 + Z.^2 - (Z * ct1 + st1 * ...
    (X * ct2 + Y * st2)).^2) / (2 * sigmaPSF^2)) .* ( ...
    (erf(abs(tmp1) / c0) .* tmp1) ./ abs(tmp1) + ...
    (erf(abs(tmp2) / c0) .* tmp2) ./ abs(tmp2));
