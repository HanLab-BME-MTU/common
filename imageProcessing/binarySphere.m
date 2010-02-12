function sph = binarySphere(radius)
%BINARYSPHERE creates a 3D spherical neighborhood/structuring element for morphological operations
% 
% sphereMat = binarySphere(radius)
% 
% This function genereates a logical matrix with values inside a sphere of
% the specified radius being true and those outside being false. This
% sphere can be used to perform morphological operations such is dilation
% and erosion on 3D binary matrices.
% 
% Input:
% 
%   radius - Positive scalar >= 1. The radius of the sphere to generate.
%    
% 
% Output:
% 
%   sphereMat - The resulting square logical matrix (neighborhood), of size
%   ~2*radius+1 (this is only exact if radius is an integer)
% 
% Hunter Elliott
% 2/2010
%

if nargin < 1 || isempty(radius) || radius < 1
    error('You must specify a radius >= 1 for the sphere!')
end

%So the radius matches the traditional 2D definition of a structuring
%element radius.
radius = radius+1;

%Get x,y,z coordinate matrices for distance-from-origin calculation
[xx,yy,zz] = meshgrid(-radius:radius,-radius:radius,-radius:radius);

%Return all points which are less than radius away from origin of
%neighborhood
sph = xx .^2 + yy .^2 + zz.^2 < radius^2;
