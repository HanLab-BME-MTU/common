function plotWTATrou(W)
% PLOTWTATROU(W) plots the coefficients from the A Trou Wavelet Transform.
% See WTATrou.m for details about that wavelet transform.
%
% Sylvain Berlemont, 2009

nScales = size(W, 3);

nCols = ceil(nScales / 2);

for i = 1:nScales
    subplot(2, nCols, i); imshow(W(:, :, i), []);
    title(['scale = ' num2str(i)]);
end