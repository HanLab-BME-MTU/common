function testKDTreeClosestPoint
% This function test the KDTreeClosestPoint.
%
% See KDTreeClosestPoint.m for details.
%
% Sebastien Besson, Oct 2011

% Generate random input and query points
X = rand(1000000,2);
C = rand(1000,2);

tic
[idx,d] = KDTreeClosestPoint(X,C);
toc

tic
D = createDistanceMatrix(X,C);
[d2,idx2] = min(D);
toc

assert(all(idx(:)==idx2(:)));
assert(all(d(:)==d2(:)))
