function m = meanShiftVector(u,p,h)
%MEANSHIFTVECTOR computes normal mean shift vector for input point
%
% m = meanShiftVector(u,x,h)
%
% Input:
%
%   u - point to calculate mean shift at. 1xM
%   x - data. NxM. This is usually points within some distance of u (support)
%   h - bandwidth of normal kernel. Scalar.
%
% Output:
%
%   m - the mean shift vector

%Hunter Elliott (but I just copied from Deepak's MeanShiftClustering
%function)
% 5/2013

d2 = sum( bsxfun(@minus,p,u) .^2 , 2);%Squared distances to u

w = exp(-.5 * d2 / h^2);%Weights from kernel

m = (w' * p) / sum(w) - u;%Shift vector. Difference here so it's one subtraction operation

