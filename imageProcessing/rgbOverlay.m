% Overlays a color mask on a grayscale input image.
%
% Inputs:
%         img     : Grayscale input image
%         overlay : Overlay mask. Nonzero elements replace values in channel 'ch'
%         ch      : channel number {1,2,3}

% Francois Aguet, June 2010.

function out = rgbOverlay(img, overlayMask, ch)

img = uint8(scaleContrast(img));
idx = overlayMask~=0;

out = repmat(img, [1 1 3]);
img(idx) = 255;
out(:,:,ch) = img;

ch = setdiff(1:3, ch);
img(idx) = 0;
out(:,:,ch(1)) = img;
out(:,:,ch(2)) = img;