function dm = distMat2(A,B,metric);
%DISTMAT computes the distance matrix of two point sets
%
% SYNOPSIS dm = distMat(M);
%
% INPUT:  M: an S x T matrix where the # of columns is the # of dimensions
%            and  the # of rows is the # of points
%
% c: 18/09/01   dT


[mA, nA] = size(A);
[mB, nB] = size(B);

if nA~=nB
    error('point set must be in same dim-space');
end;

if nargin < 3
    metric=eye(nA);
elseif length(metric)~= nA
    error('incorrect dimension of metric');
end

dm=zeros(mA,mB);

I=[1:mA]'*ones(1,mB);
I=I(:);
J=ones(mA,1)*[1:mB];
J=J(:);

Y = (A(I,:)-B(J,:))';
I = []; J = []; p = [];  % no need for I J p any more.
dm(:)=sqrt(diag(Y'*metric*Y))';