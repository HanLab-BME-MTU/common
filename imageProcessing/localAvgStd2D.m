% Computes the local average and standard deviation within a square window of side 'w'.

% Francois Aguet
% Last modified on 10/16/2009

function [avg sigma] = localAvgStd2D(image, w)

if (isEven(w))
    error('The window length w should be an odd integer.');
end;

b = (w-1)/2;
image = padarray(image, [b b], 'replicate');

h = ones(1,w);
E = conv2(h/w, h, image, 'valid');
E2 = conv2(h, h, image.^2, 'valid');

sigma = E2 - E.^2;
sigma = sqrt(sigma/(w*w - 1));
avg = E/w;