function [out] = filterGauss3D(input, sigma, borderCondition)
% filterGauss2D :	filters a data volume with a 3-D Gaussian mask
%
%     out = filterGauss3D(image, sigma, borderCondition);
%
%    INPUT: image           : 3-D input array
%           sigma           : standard deviation of the Gaussian
%           borderCondition : input for 'padarrayXT'. Default: 'symmetric'
%                             Options: 'symmetric', 'replicate', 'circular', 'antisymmetric', or a constant value
%
%    OUTPUT: out : filtered volume
%
% Francois Aguet, added 01/21/2010

if nargin < 3 || isempty(borderCondition)
    borderCondition = 'symmetric';
end

w = ceil(3*sigma); % cutoff radius of the gaussian kernel
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
g = g/sum(g);

out = convn(padarrayXT(input, [w w w], borderCondition), g', 'valid');
out = convn(out, g, 'valid');
out = convn(out, reshape(g, [1 1 2*w+1]), 'valid');