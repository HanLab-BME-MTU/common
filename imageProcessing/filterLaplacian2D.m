function [out, FL] = filterLaplacian2D(image, sigma, FL)
% filterLaplacian2D : filters an image with a 2-D Laplacian mask. Unlike
% the builtin Matlab function fspecial('Laplacian'), this function computes
% the exact laplacian convolution, using FFTs. In addition, it gives the
% ability to pass the Fourier transform of the filter for multiple calls of
% the function, saving one fft2 call per use.
%
%    out = filterLaplacian2D(image, sigma, FL);
%
%    INPUT: image           : 2-D input array
%           sigma           : standard deviation of the Gaussian
%           FL              : Fourier transform of the 2-D Laplacian mask
%                             (optional)
%
%    OUTPUT: out : filtered image
%            FL  : Fourier transform of the 2-D Laplacian mask
%
% Sylvain Berlemont, added 07/01/2010

[nrows ncols] = size(image);
hside = 5 * ceil(sigma);

if nargin < 3 || isempty(FL)
    [X Y] = meshgrid(-hside:hside);
    L = (X.^2 + Y.^2 - 2 * sigma^2) / sigma^2 .* exp(-(X.^2 + Y.^2) / (2 * sigma^2));
    L = -L / sqrt(sum(L(:).^2));
    p = nextpow2(nrows + size(L,1));
    q = nextpow2(ncols + size(L,2));

    L(2^p,2^q) = 0;
    FL = fft2(L);
else
    p = log2(size(FL,1));
    q = log2(size(FL,2));
    
    if p ~= nextpow2(nrows + 2 * hside + 1) || q ~= nextpow2(ncols + 2 * hside + 1)
        error('Inconsistent size between sigma and FL.');
    end
end

image(2^p,2^q) = 0;
FI = fft2(image);
    
out = real(ifft2(FI.*FL));
out = out(hside+1:hside+nrows, hside+1:hside+ncols);
