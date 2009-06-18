function plotTransparent(x,y,width,color,opacity)

% 
% plotTransparent(x,y,width,plotStyle)
% 
%
% Plots the input x-y data as a tansparent line with variable width. 
% 
% Inputs:
% 
%     x,y     -   The coordinates of the data to be plotted. x is optional.
%     
%     width   -   The width of the line to be plotted. Scalar or vector. If
%     scalar, the line is the same width everywhere. If a vector of same
%     length as x & y, then the width will be variable at each point.
%     Optional.
%     
%     color - The color to fill the line with. Can be a string color name
%     or a 1x3 RGB triplet. Optional.
%     
%     opacity - How transparent the line is. If 1, the line is opaque. If
%     0, the line is invisible.
% 
% Hunter Elliott, 6/2009
% 

if isempty(x)
    x = 1:length(y);
end

if nargin < 3 || isempty(width)
    width = nanstd(y) / 20;
end    

if nargin < 4 || isempty(color)
    color = 'b';
end

if nargin < 5 || isempty(opacity)
    opacity = .4;
end

n = length(x);

if length(y) ~= n
    error('X and Y must be the same length!!')
end

if length(width) == 1
    width = repmat(width,1,n);
end

if length(width) ~= n
    error('Width must be a scalar or a vector of the same length as x & y!!')
end

if any(isnan([x(:)' y(:)']))
    disp('Warning: NaN values will be ignored!')
    
    notNan = ~isnan(x(:)) & ~isnan(y(:));
    width = width(notNan);
    x = x(notNan);
    y = y(notNan);
end



xFill = [x x(end:-1:1)];
yFill = [(y - width) (y(end:-1:1) + width)];

fill(xFill,yFill,color,'EdgeColor','none','FaceAlpha',opacity)


