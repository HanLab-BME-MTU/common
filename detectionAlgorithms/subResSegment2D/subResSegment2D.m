function F = subResSegment2D(x, y, amp, sigma, l, theta, xRange, yRange, nzIdx)
% sub-resolution 2D segment model defined by 6 parameters:
%    xy      : position of the segment's center
%    amp     : mean amplitude along the segment
%    sigma   : half width of the segment
%    l       : length
%    theta   : orientation [-pi/2, pi/2)
%
% F = subResSegment2D(x, y, amp, sigma, l, theta, xRange, yRange, nzIdx)
%
% parameters:
%
% (x,y)              position of the segment's center
%
% amp                amplitude of the segment
%
% sigma              half width of the segment
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% (xRange, yRange)   2 vectors representing the 2-dimensional support of
%                    the segment. This support can be determined using
%                    subResSegment2DSupport() function.
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

xRange = xRange - x;
yRange = yRange - y;

N = numel(yRange);
M = numel(xRange);

if nargin < 10 || isempty(nzIdx)
    nzIdx = 1:N*M;
end

[X Y] = meshgrid(xRange, yRange);
X = X(nzIdx);
Y = Y(nzIdx);

ct = cos(theta);
st = sin(theta);

C0 = (1/2).*amp.*erf(2.^(-1/2).*l.*sigma.^(-1)).^(-1);
C1 = (1/2).*2.^(-1/2).*sigma.^(-1);

F = zeros(N,M);

F(nzIdx) = C0 * exp((-1/2).*sigma.^(-2).*(Y.*ct-X.*st).^2).*(...
    erf(C1.*(l+2.*X.*ct+2.*Y.*st))+...
    erf(C1.*(l-2*X*ct-2*Y*st)));
