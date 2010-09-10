function dF = subResSegment2DJacobian(x, y, amp, sigma, l, theta, xRange, yRange, nzIdx, paramSelector)
% Partial derivatives of a sub-resolution 2D segment against its parameters.
% dF = subResSegment2DJacobian(x, y, amp, sigma, l, theta, xRange, yRange, nzIdx, paramSelector)
%
% parameters:
% (x,y)              position of the segment's center
%
% amp                amplitude of the segment
%
% sigma              half width of the gaussian PSF model.
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
%                    M = numel(xRange)) where the segment is not null. If
%                    not provided, nzIdx = ones(N*M, 1). These indices can
%                    be determined using subResSegment2DSupport() function.
%
% paramSelector      
%
% output:
% dF is a NxMxP matrix.
%
% Sylvain Berlemont, 2010

xRange = xRange - x;
yRange = yRange - y;

N = numel(yRange);
M = numel(xRange);

if nargin < 7 || isempty(nzIdx)
    nzIdx = ones(N*M, 1);
end

if nargin < 8 || isempty(paramSelector)
    paramSelector = ones(7,1);
end
    
P = nnz(paramSelector);

ct = cos(theta);
st = sin(theta);

[X Y] = meshgrid(xRange, yRange);
X = X(nzIdx);
Y = Y(nzIdx);

s_1 = sigma^(-1);
s_2 = sigma^(-2);
s_3 = sigma^(-3);

C1 = (1/2).*exp((-1/2)*s_2.*(Y*ct-X*st).^2);
C2 = erf(l / (2^(1/2) * sigma));
C3 = erf((1/2).*2.^(-1/2)*s_1.*(l+2*X*ct+2*Y*st));
C4 = erf((1/2).*2.^(-1/2)*s_1.*(l-2*X*ct-2*Y*st));
C5 = exp((-1/8)*s_2.*(l-2*X*ct-2*Y*st).^2);
C6 = exp((-1/8)*s_2.*(l+2*X*ct+2*Y*st).^2);
C7 = exp((-1/2).*l.^2*s_2);
C8 = (2.*pi.^(-1)).^(1/2);
C9 = (2.*pi).^(-1/2);

dFuncs = {@dFdXc, @dFdYc, @dFdA, @dFds, @dFdl, @dFdt};

dF = zeros(N,M,P);

paramSet = 1:numel(paramSelector);

offset = 0;

for iParam = paramSet(paramSelector)
    dF(nzIdx + offset) = dFuncs{iParam}();
    
    offset = offset + N * M;
end

    function res = dFdXc
        res = amp*C1*s_2.*C2.^(-1).*(C5*C8*sigma*ct-C6*C8*sigma*ct+...
            st*(C3+C4).*(X*st-Y*ct));
    end

    function res = dFdYc
        res = amp*C1*s_2.*C2.^(-1).*(C5*C8*sigma*st-C6*C8*sigma*st+...
            ct*(C3+C4).*(Y*ct-X*st));
    end

    function res = dFdA
        res = C1.*C2.^(-1).*(C3+C4);
    end

    function res = dFds
        res = amp*C1.*C2.^(-2).*(C7.*l.*C8*s_2.*(C3+C4)+s_3.*C2.*(C3+C4).*...
            (Y*ct-X*st).^2+C9*s_2.*C2.*(C5.*(-l+2*X*ct+2*Y*st)-C6.*...
            (l+2*X*ct+2*Y*st)));
    end

    function res = dFdl
        res = amp*C1.*C9*s_1.*C2.^(-2).*((C6+C5).*C2-2*C7.*(C3+C4));
    end

    function res = dFdt
        res = amp*C1.*C2.^(-1).*(C5.*C9.*s_1.*(2*X*st-2*Y*ct)+C6.*C8*s_1.*...
            (Y*ct-X*st)+s_2*(C3+C4).*(Y*ct-X*st).*(X*ct+Y*st));
    end

%     function res = dFdBg
%         res = ones(N,M);
%     end
end