% Overlays a color mask on a grayscale input image.
%
% Inputs:
%         img     : Grayscale input image
%         overlay : Overlay mask. Nonzero elements replace values in channel 'ch'
%         ch      : channel number {1,2,3}

% Francois Aguet, June 2010.

function img = rgbOverlay(img, overlay, ch)

img = repmat(uint8(scaleContrast(img)), [1 1 3]);
c = img(:,:,ch);
c(overlay==max(overlay(:))) = 255;
img(:,:,ch) = c;