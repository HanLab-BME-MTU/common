function h = plotTransparent(varargin)

%
% h = plotTransparent(x,y,width,color,opacity,perpendicularOffset,useLimits)
% h = plotTransparent(ah,...
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
%     correctAspectRatio - if 1, expected axes limits are used to adjust the width
%                 of the line in x/y so that it doesn't appear deformed
%                 without axes equal. If an axes handle has been given, the
%                 dataAspectRatio of these axes will be used to calculate
%                 the correction.
%                 Default: false.
%                 THIS OPTION DOES NOT WORK PROPERLY AT THE MOMENT
%
%     ah : handles of axes into which to plot
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
%           x = [1,2,3,4,4,6];y=[1,2,1,1,2,1];
%           plotTransparent(x,y,0.1) 
%           axis equal % if this is omitted, the plot will not look well.
%
% Hunter Elliott, 6/2009
% adapted for any 2D lines by Jonas


if nargin  == 0
    error('plotTransparent needs at least one nonempty input argument')
end

if isscalar(varargin{1}) && ishandle(varargin{1}) &&  strcmp(get(varargin{1},'Type'),'axes')
    ah = varargin{1};
    varargin(1) = [];
else
    ah = [];
end

% overwrite nargin
numArgIn = length(varargin);

if numArgIn  == 0
    error('plotTransparent needs at least one nonempty input argument')
end

x = varargin{1}(:);

if numArgIn < 2 || isempty(varargin{2})
    if isempty(x)
        error('plotTransparent needs at least one nonempty x or y')
    end
    y = x;
    x = [];
else
    y = varargin{2}(:);
end

if isempty(x)
    x = (1:length(y))';
end

if numArgIn < 3 || isempty(varargin{3})
    width = nanstd(y) / 20;
else
    width = varargin{3}(:);
end

if numArgIn < 4 || isempty(varargin{4})
    color = 'b';
else
    color = varargin{4};
end

if numArgIn < 5 || isempty(varargin{5})
    opacity = .4;
else
    opacity = varargin{5};
end

if numArgIn < 6 || isempty(varargin{6})
    perpendicularOffset = true;
else
    perpendicularOffset = varargin{6};
end

if numArgIn < 7 || isempty(varargin{7})
    useLimits = false;
else
    useLimits = varargin{7};
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
    % Hunter's version - updated by Jonas
    offset = [zeros(size(x,1),1),width];
    
else
    % Jonas' version
    
    xy = [x,y];
    
    
   
    
    % find derivative to have proper corners
    [n_xy,e_xy] = normList(diff(xy));
    
     
    
    
    
    offset = [-e_xy(:,2),e_xy(:,1)];
    
    % check whether to use limits
    if useLimits
        if isempty(ah)
            fh = figure('visible','off');
            tmpAx = axes('Parent',fh);
            plot(x,y);
            aspectRatio = get(tmpAx,'DataAspectRatio');
            delete(fh)
        else
            aspectRatio = get(ah,'DataAspectRatio');
        end
        stretchRatio = aspectRatio(1)/aspectRatio(2);
        % offset gets stretched in x vs y. However, the line gets stretched
        % in x vs y as well. Therefore, in order to look perpendicular,
        % both x and y have to be corrected
        offset(:,1) = offset(:,1) * (stretchRatio);
        offset(:,2) = offset(:,2) / (stretchRatio);
        % this gives the right angles, but unfortunately, the correction is
        % wrong. However, I am unable to figure this out right now.
    end
    
    % offset should be average of normal to the line segments at
    % intersections. In order to get constant line thickness, the averaged
    % offset vector needs to be stretched by 1/sin(alpha), otherwise, there
    % is an apparent narrowing of the line, especially at acute angles
    tmp = [offset(1,:);(offset(2:end,:) + offset(1:end-1,:))/2;offset(end,:)];
    tmp(2:end-1,:) = tmp(2:end-1,:) ./...
        (dot(offset(1:end-1,:),tmp(2:end-1,:),2) * [1,1]);
    % Adjust with based on user input. Not using bsxfun here because of LCCB
    offset = (width * [1 1]) .* tmp;
    
    
    
end
    
    
    
    % create x,y for fill (could also use multiple patches if multiple colors
    % are needed 
    xFill = [x + offset(:,1); x(end:-1:1) - offset(end:-1:1,1)];
    yFill = [y + offset(:,2); y(end:-1:1) - offset(end:-1:1,2)];

    if ~isempty(ah)
        set(gca,'DataAspectRatio',aspectRatio,'NextPlot','add')
    end
ph = fill(xFill,yFill,color,'EdgeColor','none','FaceAlpha',opacity);
% hold on, plot(xFill,yFill,'k')
% plot(x,y,'ok')

if nargout > 0
    h = ph;
end