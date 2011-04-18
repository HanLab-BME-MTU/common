% Overlays a color mask on a grayscale input image.
%
% Inputs:
%         img          : Grayscale input image.
%         mask         : Overlay mask. Nonzero elements are colored.
%         overlayColor : 3-element RGB vector.
%         {iRange}     : dynamic range of 'img'.

% Francois Aguet, April 2011.

function imgRGB = rgbOverlay(img, mask, overlayColor, iRange)

if nargin<4
    iRange = [];
end

[chR chG chB] = deal(scaleContrast(img, iRange));
maskIdx = mask~=0;
chR(maskIdx) = chR(maskIdx)*overlayColor(1);
chG(maskIdx) = chG(maskIdx)*overlayColor(2);
chB(maskIdx) = chB(maskIdx)*overlayColor(3);
imgRGB = uint8(cat(3, chR, chG, chB));