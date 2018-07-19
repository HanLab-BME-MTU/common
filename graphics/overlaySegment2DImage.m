function overlaySegment2DImage(ima, segmentParams)
% Overlay a set of lines on an image, where each line represents a
% sub-resolution 2D segment.
%
% overlaySegment2DImage(ima, segmentParams)
%
% INPUT:
% ima             : the image
%
% segmentParams   : nx5 matrix where n is the number of segments and their
%                   parameters, i.e. xC, yC, A, l, t are stored column-wise.
%
% Sylvain Berlemont, 2010

imagesc(ima),colormap gray,axis image,axis off;

hold on;

xC = segmentParams(:,1);
yC = segmentParams(:,2);
l = segmentParams(:,5);
t = segmentParams(:,6);

ct = cos(t);
st = sin(t);

line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', ...
     [yC - (l / 2) .* st, yC + (l / 2) .* st]', ...
    'Color', 'g');

line(xC, yC, 'Color', 'g', 'Marker', '.', 'LineStyle', 'none');
