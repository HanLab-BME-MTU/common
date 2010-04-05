function [dFdXc dFdYc dFdA dFds dFdl dFdt] = ...
    dLSegment2DJacobian(xRange, yRange,A, sigmaPSF, l, theta, nzIdx)
% Partial derivatives of a diffraction-limited 2D segment against each of
% its parameter.
% [dFdXc dFdXy dFdA dFds dFdl dFdt] =
%           dLSegment2DJacobian(xRange, yRange, A, sigmaPSF, l, theta)
%
% parameters:
% (xRange, yRange)   2 vectors representing the 2-dimensional support of
%                    the segment. This support can be determined using
%                    dLSegment2DSupport() function.
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
%                    M = numel(xRange)) where the segment is not null. If
%                    not provided, nzIdx = ones(N*M, 1). These indices can
%                    be determined using dLSegment2DSupport() function.
%
% output:
% dFdXc dFdXy dFdA dFds dFdl dFdt are NxM matrices. 
%
% Sylvain Berlemont, 2010

N = numel(yRange);
M = numel(xRange);

if nargin < 7 || isempty(nzIdx)
    nzIdx = ones(N*M, 1);
end

ct = cos(theta);
st = sin(theta);

[X Y] = meshgrid(xRange, yRange);
X = X(nzIdx);
Y = Y(nzIdx);

s_1 = sigmaPSF^(-1);
s_2 = sigmaPSF^(-2);
s_3 = sigmaPSF^(-3);

C1 = (1/2).*exp((-1/2)*s_2.*(Y*ct-X*st).^2);
C2 = erf(l / (2^(1/2) * sigmaPSF));
C3 = erf((1/2).*2.^(-1/2)*s_1.*(l+2*X*ct+2*Y*st));
C4 = erf((1/2).*2.^(-1/2)*s_1.*(l-2*X*ct-2*Y*st));
C5 = exp((-1/8)*s_2.*(l-2*X*ct-2*Y*st).^2);
C6 = exp((-1/8)*s_2.*(l+2*X*ct+2*Y*st).^2);
C7 = exp((-1/2).*l.^2*s_2);
C8 = (2.*pi.^(-1)).^(1/2);
C9 = (2.*pi).^(-1/2);

dFdXc = zeros(N,M);
dFdXc(nzIdx) = A*C1*s_2.*C2.^(-1).*(C5*C8*sigmaPSF*ct-C6*C8*sigmaPSF*ct+...
    st*(C3+C4).*(X*st-Y*ct));

dFdYc = zeros(N,M);
dFdYc(nzIdx) = A*C1*s_2.*C2.^(-1).*(C5*C8*sigmaPSF*st-C6*C8*sigmaPSF*st+...
    ct*(C3+C4).*(Y*ct-X*st));

dFdA = zeros(N,M);
dFdA(nzIdx) = C1.*C2.^(-1).*(C3+C4);

dFds = zeros(N,M);
dFds(nzIdx) = A*C1.*C2.^(-2).*(C7.*l.*C8*s_2.*(C3+C4)+s_3.*C2.*(C3+C4).*...
    (Y*ct-X*st).^2+C9*s_2.*C2.*(C5.*(-l+2*X*ct+2*Y*st)-C6.*(l+2*X*ct+2*Y*st)));

dFdl = zeros(N,M);
dFdl(nzIdx) = A*C1.*C9*s_1.*C2.^(-2).*((C6+C5).*C2-2*C7.*(C3+C4));

dFdt = zeros(N,M);
dFdt(nzIdx) = A*C1.*C2.^(-1).*(C5.*C9.*s_1.*(2*X*st-2*Y*ct)+C6.*C8*s_1.*...
    (Y*ct-X*st)+s_2*(C3+C4).*(Y*ct-X*st).*(X*ct+Y*st));
