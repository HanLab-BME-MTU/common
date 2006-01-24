function dm = distMat2(A,B,metric);
%DISTMAT2 computes the distance matrix of two point sets
%
% SYNOPSIS dm = distMat2(A,B,metric);
%
% INPUT:  A, B: Arrays containg rowA, rowB points with nCol dimensions.
%         metric (optional): weight of the dimensions
%
% OUTPUT: dm: rowA-by-rowB array of distances between points in A and B
%
% c: 18/09/01   dT


[mA, nA] = size(A);
[mB, nB] = size(B);

if nA~=nB
    error('point set must be in same dim-space');
end;

if nargin < 3
    metric=ones(1,nA);
    metricIsSquare = 0;
elseif ~all(size(metric) == [nA,nB])
    error('incorrect dimension of metric');
else
    metricIsSquare = 1;
end

dm=zeros(mA,mB);

[I,J]=ndgrid(1:mA,1:mB);
I=I(:);
J=J(:);

% repeat the entries in A and B so that we'll take the difference between
% all of them
Y = (A(I,:)-B(J,:))';

% calculate distance as the rood of the squared sum. If metric is not a
% square, we can use a memory-saving calculation
if metricIsSquare
    dm(:)=sqrt(diag(Y'*metric*Y))';
else
    dm(:) = sqrt(sum(Y'.^2.*repmat(metric,mA*mB,1),2));
end