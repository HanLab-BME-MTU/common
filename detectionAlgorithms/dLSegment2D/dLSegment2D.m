function F = subResSegment2D(xRange, yRange, A, sigmaPSF, l, theta, nzIdx)
% sub-resolution 2D segment model.
% F = subResSegment2D(xRange, yRange, A, sigmaPSF, l, theta, nzIdx)
%
% parameters:
% (xRange, yRange)   2 vectors representing the 2-dimensional support of
%                    the segment. This support can be determined using
%                    subResSegment2DSupport() function.
%
% A                  amplitude of the segment
%
% sigmaPSF           half width of the gaussian PSF model.
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% nzIdx              linear indices of a NxM matrix (N = numel(yRange) and
%                    M = numel(xRange)) where the model is defined. If not
%                    provided, nzIdx = 1:N*M. These indices can be
%                    determined using subResSegment2DSupport() function.
%
% output:
% F                  the model defined on a NxM matrix.
%
% Sylvain Berlemont, 2010

N = numel(yRange);
M = numel(xRange);

if nargin < 7 || isempty(nzIdx)
    nzIdx = 1:N*M;
end

[X Y] = meshgrid(xRange, yRange);
X = X(nzIdx);
Y = Y(nzIdx);

ct = cos(theta);
st = sin(theta);

C0 = (1/2).*A.*erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-1);
C1 = (1/2).*2.^(-1/2).*sigmaPSF.^(-1);

F = zeros(N,M);

F(nzIdx) = C0 * exp((-1/2).*sigmaPSF.^(-2).*(Y.*ct-X.*st).^2).*(...
    erf(C1.*(l+2.*X.*ct+2.*Y.*st))+...
    erf(C1.*(l-2*X*ct-2*Y*st)));
