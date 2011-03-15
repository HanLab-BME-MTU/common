%plotScaleBar(width, h, varargin) adds a scale bar to a figure.
%
% Inputs:  width : width of the scale bar, in x-axis units.
%              h : axes handle. If empty, current axes ('gca') are used.
%       varargin : optional inputs, always in name/value pairs:
%                  'Location' : {'NorthEast', 'SouthEast', 'SouthWest', 'NorthWest'}
%                  'Label' : string
%                  'FontName'
%                  'FontSize'
%
% Example: plotScalebar(500, [], 'Label', '500 nm', 'Location', 'SouthEast');


% Note: there is a bug in the 'patch' function of Matlab2010b, hence the 0.999 factor.
%
% Francois Aguet, March 14 2011


function plotScaleBar(width, h, varargin)

% fixed parameters
b = 40;


if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

if nargin<2 || isempty(h)
    h = gca;
end


idx = find(strcmpi(varargin, 'location'));
if ~isempty(idx)
    location = lower(varargin{idx+1});
else
    location = 'southwest';
end

idx = find(strcmpi(varargin, 'Label'));
if ~isempty(idx)
    label = varargin{idx+1};
else
    label = [];
end

idx = find(strcmpi(varargin, 'FontName'));
if ~isempty(idx)
    fontName = varargin{idx+1};
else
    fontName = 'Helvetica';
end

idx = find(strcmpi(varargin, 'FontSize'));
if ~isempty(idx)
    fontSize = varargin{idx+1};
else
    fontSize = 12;
end


XLim = get(h, 'XLim');
YLim = get(h, 'YLim');

lx = diff(XLim);
ly = diff(YLim);

height = ly/b; % height of the scale bar
dx = lx/b; % offset from image border

if ~isempty(label)
    % get height of default text bounding box
    h = text(0, 0, ' ', 'FontName', fontName, 'FontSize', fontSize);
    extent = get(h, 'extent');
    textHeight = extent(4);
    delete(h);
else
    textHeight = 0;
end

hold on;
set(gcf, 'InvertHardcopy', 'off');
switch location
    case 'northeast'
        fill([lx-width lx lx lx-width]-dx, [height height 0 0]+dx,...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(lx-dx-width/2, dx+height, label, 'Color', 'w',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Center',...
                'FontName', fontName, 'FontSize', fontSize);
        end
    case 'southeast'
        patch([lx-width lx lx lx-width]-dx, [ly ly ly-height ly-height]-max(dx,textHeight),...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(lx-dx-width/2, ly-max(dx,0.9*textHeight), label, 'Color', 'w',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Center',...
                'FontName', fontName, 'FontSize', fontSize);
        end
    case 'southwest'
        fill([0 width width 0]+dx, [ly ly ly-height ly-height]-max(dx,textHeight),...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(dx+width/2, ly-max(dx,0.9*textHeight), label, 'Color', 'w',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Center',...
                'FontName', fontName, 'FontSize', fontSize);
        end
    case 'northwest'
        fill([0 width width 0]+dx, [height height 0 0]+dx,...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(dx+width/2, dx+height, label, 'Color', 'w',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Center',...
                'FontName', fontName, 'FontSize', fontSize);
        end
end
