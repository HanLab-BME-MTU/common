function dLSegment2DPlot(I, params)
% dLSegment2DPlot(I, params)
% parameters:
%
% I               image
%
% params          nx6 matrix where n is the number of segments and their
%                 parameters, i.e. xC, yC, A, l, t are stored column-wise.
%

imshow(I, []); hold on;

xC = params(:,1);
yC = params(:,2);
l = params(:,4);
t = params(:,5);

line([xC, xC + (l / 2) .* cos(t)], [yC, yC + (l / 2) .* sin(t)],'r');
line([xC, xC + (l / 2) .* cos(t + pi)], [yC, yC + (l / 2) .* sin(t + pi)], 'r');

plot(xC, yC, 'r.');