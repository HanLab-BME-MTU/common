% filterLoG : filters an image with a 2-D Laplacian of Gaussian filter. Unlike
% the built-in Matlab function fspecial('Laplacian'), this function computes
% the exact convolution, using FFTs.
%
%    out = filterLoG(img, sigma);
%
%    INPUT: img   : 2-D input image
%           sigma : standard deviation of the Gaussian
%
%    OUTPUT: y : LoG filtered image
%
% Francois Aguet, June 28, 2010

function y = filterLoG(img, sigma)

[ny,nx] = size(img);

[w1,w2] = meshgrid(-nx/2:nx/2-1, -ny/2:ny/2-1);
w1 = fftshift(w1);
w2 = fftshift(w2);

w1 = w1*2*pi/nx;
w2 = w2*2*pi/ny;
I = fft2(img);
G = 4*pi*sigma^2 * (w1.^2 + w2.^2) .* exp(-0.5*sigma^2*(w1.^2 + w2.^2));
y = real(ifft2(I.*G));