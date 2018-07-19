function [ X, Y, Z, pxInd, varargout] = randmask( mask, n, m, varargin)
%RANDMASK Find random points located in a ND pixel/voxel mask on a
%continuous basis
%
% INPUT
% mask - logical (binary) 2D or 3D mask
% n - number of points
%
% OUTPUT
% X - 1st dimension (column)
% Y - 2nd dimension (row)
% Z - 3rd dimension relative to mask
% pxInd - linear index of the pixel chosen in the mask
% varargout - 4th dimension and beyond
%
% ALGORITHM
%
% The function selects a random integer to first select a pixel within a
% pixel based mask. Then with each pixel selected a random continuous value
% is selected.
%
% This is an alternative to selecting more random values than needed and
% rejecting ones outside of the mask.
%
% The strength of this approach is that this should work for an arbitrary
% pixel defined mask including in n-dimensions.
%
% Example
%
% [X,Y] = meshgrid(1:101,1:101);
% mask = zeros(size(X));
% mask(hypot(X-50,Y-50) > 30) = 1;
% [rmX,rmY] = randmask(mask,1e3,1);
% figure;
% imshow(mask,[]);
% hold on;
% plot(rmX,rmY,'.');
%
% See also randcirc, inpolygon, unifrnd


% Mark Kittisopikul, Ph.D.
% Northwestern
% October 2017

if(nargin < 2)
    n = 1;
end

if(nargin < 3)
    m = n;
end

ndim = max(nargout,2);
if(nargout > 3)
    % The 4th output is the pixel index,
    % so reduce the number of dimensions by one
    ndim = ndim - 1;
end

dims = [n m varargin];
colons = {':'};
colons = colons(ones(size(dims)));

% rand is in the open interval (0,1)
% compensate by making this (0,1] ??
r = rand(dims{:},ndim)./(1-eps(1));

ind = find(mask);

% First randomly select a pixel in the mask
pxInd = randi(length(ind),dims{:});

% Get the pixel coordinates
if(nargout <= 2)
    % No Z processing needed
    [Y,X] = ind2sub(size(mask),ind(pxInd));
elseif(nargout == 3 || nargout == 4)
    [Y,X,Z] = ind2sub(size(mask),ind(pxInd));
    Z = Z - r(colons{:},3) +0.5;
else
    % ND processing
    [Y,X,Z,varargout{:}] = ind2sub(size(mask),ind(pxInd));
    Z = Z - r(colons{:},3) +0.5;
    for d=5:nargout
        varargout{d-4} = varargout{d-4}  - r(:,d-2) +0.5;
    end
end

% Always handle the first two dimensions
Y = reshape(Y,dims{:}) - r(colons{:},1) +0.5;
X = reshape(X,dims{:}) - r(colons{:},2) +0.5;

end
