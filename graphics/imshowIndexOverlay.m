function imHan = imshowIndexOverlay(L,alpha,cMap)
%IMSHOWINDDEXOVERLAY shows a transparent RGB overlay of the input index image
%
% imshowIndexOverlay(L)
% imHan = imshowIndexOverlay(L,alpha,cMap)
%
% Hunter Elliott
% 4/16/2016
%

if nargin < 3 || isempty(cMap)
    cMap = randomColormap(numel(unique(L)),42);
end

if nargin < 2 || isempty(alpha)
    alpha = double(L > 0) * .5;
elseif isscalar(alpha)
    alpha = ones(size(L)) * alpha;
end

overlayIm = ind2rgb(L,cMap);

hold on
imHan = imshow(overlayIm);
imHan.AlphaData = alpha;


