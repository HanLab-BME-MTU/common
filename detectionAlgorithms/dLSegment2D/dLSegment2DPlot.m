function dLSegment2DPlot(I, params, color)
% dLSegment2DPlot(I, params)
%
% Overlay a set of lines on an image, where each line represents a
% diffraction-limited 2D segment.
%
% parameters:
%
% I               the image
%
% params          nx6 matrix where n is the number of segments and their
%                 parameters, i.e. xC, yC, A, l, t are stored column-wise.
%

%imshow(I, []); hold on;

if nargin < 3 || isempty(color) || ~ischar(color)
    color = 'g';
end

xC = params(:,1);
yC = params(:,2);
l = params(:,4);
t = params(:,5);

line([xC - (l / 2) .* cos(t), xC + (l / 2) .* cos(t)]', ...
     [yC - (l / 2) .* sin(t), yC + (l / 2) .* sin(t)]', ...
    'Color', color);

plot(xC, yC, 'Color', color, 'LineStyle', '.');
