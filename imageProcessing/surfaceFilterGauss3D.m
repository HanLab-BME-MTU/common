function [d2X,d2Y,d2Z] = surfaceFilterGauss3D(input, sigma, borderCondition)
%SURFACEFILTERGAUSS3D :	surface filters a data volume with a 3-D Gaussian second derivative kernel
%
%    [d2X,d2Y,d2Z] = surfaceFilterGauss3D(image, sigma, borderCondition);
%
%       Filters the input matrix using partial second derivatives of a gaussian,
%       giving a filtered "surface" image.
%
%    INPUT: image           : 3-D input array
%           sigma           : standard deviation of the Gaussian to use
%                             derivatives of for filtering
%           borderCondition : input for 'padarrayXT'. Default: 'symmetric'
%                             Options: 'symmetric', 'replicate', 'circular', 'antisymmetric', or a constant value
%
%    OUTPUT: [d2X,d2Y,d2Z] : Matrices filtered with partial derivatives of the
%                         gaussian in the X, Y and Z directions
%                         respectively, corresponding to matrix dimensions
%                         2, 1 and 3 respectively.
%
% Hunter Elliott, added 01/21/2010
% Modelled after filterGauss3D.m - thanks Francois!

if nargin < 3 || isempty(borderCondition)
    borderCondition = 'symmetric';
end

w = ceil(3*sigma); % cutoff radius of the gaussian kernel
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
d2g = (x.^2-1) / sigma^2 .* exp(-x.^2/(2*sigma^2));
gSum = sum(g);
g = g/gSum;
d2g = d2g/gSum;
d2g = d2g - sum(d2g);

d2X = convn(padarrayXT(input, [w w w], borderCondition), d2g, 'valid');
d2X = convn(d2X, g', 'valid');
d2X = convn(d2X,reshape(g,[1 1 2*w+1]),'valid');

d2Y = convn(padarrayXT(input, [w w w], borderCondition), g, 'valid');
d2Y = convn(d2Y, d2g', 'valid');
d2Y = convn(d2Y,reshape(g,[1 1 2*w+1]),'valid');

d2Z = convn(padarrayXT(input, [w w w], borderCondition), g, 'valid');
d2Z = convn(d2Z, g', 'valid');
d2Z = convn(d2Z,reshape(d2g,[1 1 2*w+1]),'valid');
