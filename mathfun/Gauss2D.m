function [out, G] = Gauss2D(image, sigma, borderCondition)
% Gauss2D :	filters an image with a 2-D Gaussian mask
%
%    [out, G] = Gauss2D(image, sigma, borderCondition);
%
%    INPUT: image           : 2-D input array
%           sigma           : standard deviation of the Gaussian
%           borderCondition : input for 'padarrayXT'. Default: 'symmetric'
%
%    OUTPUT: out : filtered image
%            G   : Gaussian mask
%
% Francois Aguet, added 01/21/2010

if nargin < 3 || isempty(borderCondition)
    borderCondition = 'symmetric';
end

w = ceil(3*sigma); % cutoff radius of the gaussian kernel
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
G = g'*g;

out = conv2(g', g, padarrayXT(image, [w w], borderCondition), 'valid');