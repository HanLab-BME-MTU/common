function h = plotTransparent(x,y,width,color,opacity,perpendicularOffset)

%
% h = plotTransparent(x,y,width,plotStyle)
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
%     perpendicularOffset - if 1, the code assumes equal axes and thus calculates the
%                 line such that it is expanded at a right angle to the
%                 individual line segments. This only looks good if you set
%                 'axis equal' or if you plot onto an image created with
%                 imshow, but it allows to plot any 2D line, not just those
%                 that can be described as y = f(x). Optional. Default: 1
%
%                 
%
% Output: h   -    Handle to the patch object
%
% Example 
%
%           phi = -pi:0.01:pi;
%           figure,
%           plotTransparent(sin(phi),cos(phi)),axis equal
%
% Hunter Elliott, 6/2009
% adapted for any 2D lines by Jonas

% note: it turns out that Hunter's version looks equally crappy with
% unequal axes (example: figure,plotTransparent(rand(25,1),[],[],[],[],0) )
% Therefore, I set the perpendicularOffset as default.


if nargin  == 0
    error('plotTransparent needs at least one nonempty input argument')
end

if nargin < 2 || isempty(y)
    if isempty(x)
        error('plotTransparent needs at least one nonempty x or y')
    end
    y = x;
    x = [];
end

if isempty(x)
    x = 1:length(y);
end

x = x(:);
y = y(:);

if nargin < 3 || isempty(width)
    width = nanstd(y) / 20;
end

if nargin < 4 || isempty(color)
    color = 'b';
end

if nargin < 5 || isempty(opacity)
    opacity = .4;
end

if nargin < 6 || isempty(perpendicularOffset)
    perpendicularOffset = true;
end

n = length(x);

if length(y) ~= n
    error('X and Y must be the same length!!')
end

if length(width) == 1
    width = repmat(width,n,1);
end

width = width(:);

if length(width) ~= n
    error('Width must be a scalar or a vector of the same length as x & y!!')
end

if any(isnan([x; y]))
    disp('Warning: NaN values will be ignored!')
    
    notNan = ~isnan(x) & ~isnan(y);
    width = width(notNan);
    x = x(notNan);
    y = y(notNan);
end



if ~perpendicularOffset
    % Hunter's version
    xFill = [x(:)' x(end:-1:1)'];
    yFill = [(y(:)' - width(:)') (y(end:-1:1)' + width(end:-1:1)') ];
    
else
    % Jonas' version
    
    xy = [x,y];
    % find derivative to have proper corners
    [n_xy,e_xy] = normList(diff(xy));
    offset = [-e_xy(:,2),e_xy(:,1)];
    % offset should be average of normal to the line segments at intersections
    % not using bsxfun here because of LCCB
    offset = (width * [1 1]) .* [offset(1,:);(offset(2:end,:) + offset(1:end-1,:))/2;offset(end,:)];
    % create x,y for fill (could also use multiple patches if multiple colors
    % are needed
    
    xFill = [xy(:,1) + offset(:,1); xy(end:-1:1,1) - offset(end:-1:1,1)];
    yFill = [xy(:,2) + offset(:,2); xy(end:-1:1,2) - offset(end:-1:1,2)];
end

ph = fill(xFill,yFill,color,'EdgeColor','none','FaceAlpha',opacity);

if nargout > 0
    h = ph;
end