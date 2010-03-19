function F = dLSegment2D(xRange, yRange, xC, yC, A, Bg, sigmaPSF, l, theta)
% 2D Diffraction-limited Segment Model
% F = dLSegment2D(xRange, yRange, xC, xC, A, Bg, sigmaPSF, l, theta)
%
% parameters:
% (xRange, yRange)   2 vectors representing the 2-dimensional domain (e.g.
%                    xRange = -10:.1:10, yRange = -5:.1:5
%
% (xC,yC)            center of the segment
%
% A                  amplitude of the segment
%
% Bg                 amplitude of the background
%
% sigmaPSF           half width of the gaussian PSF model.
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% output:
% F is a NxM matrix where N = numel(X) and M = numel(Y).
%
% Sylvain Berlemont, 2009

ct = cos(theta);
st = sin(theta);

[X Y] = meshgrid(xRange, yRange);

X = X - xC;
Y = Y - yC;

C0 = (1/2).*A.*erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-1);
C1 = (1/2).*2.^(-1/2).*sigmaPSF.^(-1);

F = Bg + C0 * exp((-1/2).*sigmaPSF.^(-2).*(Y.*ct-X.*st).^2).*(...
    erf(C1.*(l+2.*X.*ct+2.*Y.*st))+...
    erf(C1.*(l-2*X*ct-2*Y*st)));
