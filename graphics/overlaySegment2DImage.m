function overlaySegment2DImage(I, segmentParams)
% Overlay a set of lines on an image, where each line represents a
% sub-resolution 2D segment.
%
% overlaySegment2DImage(I, segmentParams)
%
% INPUT:
% I               : the image
%
% segmentParams   : nx5 matrix where n is the number of segments and their
%                   parameters, i.e. xC, yC, A, l, t are stored column-wise.
%
% Sylvain Berlemont, 2010

imagesc(I),colormap gray,axis image,axis off;

hold on;

xC = segmentParams(:,1);
yC = segmentParams(:,2);
l = segmentParams(:,4);
t = segmentParams(:,5);

ct = cos(t);
st = sin(t);

line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', ...
     [yC - (l / 2) .* st, yC + (l / 2) .* st]', ...
    'Color', 'g');

line(xC, yC, 'Color', 'g', 'LineStyle', '.');
