function awtDisplay(W)
% awtDisplay(W) plots the coefficients from the A Trou Wavelet Transform.
% See WTATrou.m for details about that wavelet transform.
%
% Sylvain Berlemont, 2009

n = size(W, 3);

nCols = ceil(n / 2);

for i = 1:n-1
    subplot(2, nCols, i); imshow(W(:, :, i), []);
    title(sprintf('W_{%d}', i));
end

subplot(2, nCols, n); imshow(W(:, :, n), []);
title(sprintf('A_{%d}', n-1));