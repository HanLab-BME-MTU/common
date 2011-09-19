% Computes the local average and standard deviation within a square window of side 'w'.

% Francois Aguet
% Last modified on 09/19/2011

function [avg sigma] = localAvgStd2D(img, w)

if mod(w+1, 2)
    error('The window length w should be an odd integer.');
end;

nanMask = isnan(img);

b = (w-1)/2;

% kernel
h = ones(1,w);

% count of non-NaN elements
n = conv2(h, h, padarray(double(~nanMask), [b b], 'replicate'), 'valid');
img(nanMask) = 0;

img = padarray(img, [b b], 'replicate');
E = conv2(h, h, img, 'valid');
E2 = conv2(h, h, img.^2, 'valid');

sigma = E2 - E.^2./n;
sigma(sigma<0) = 0;
sigma = sqrt(sigma./(n - 1));
avg = E./n;

avg(nanMask) = NaN;
sigma(nanMask) = NaN;