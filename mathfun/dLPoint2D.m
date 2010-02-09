function F = dLPoint2D(X, Y, xc, yc, A, sigmaPSF)
% Diffraction-limited Point Model
% F = dLPoint2D(X, Y, xc, xy, A, sigmaPSF);
%
% parameters:
% (X, Y)       2 vectors representing the 2-dimensional domain (e.g X =
%              -10:.1:10, Y = -5:.1:5
%
% (xc,yc)      center of the segment
%
% A            amplitude of the segment
%
% sigmaPSF     half width of the gaussian PSF model.
%
% output:
% F is a NxM matrix where N = numel(X) and M = numel(Y).
%
% Sylvain Berlemont, 2009

[X Y] = meshgrid(X,Y);
X = X - xc;
Y = Y - yc;
F = A * exp(-0.5 * (X.^2 + Y.^2) / sigmaPSF^2);
end