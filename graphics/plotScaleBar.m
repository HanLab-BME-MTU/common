%plotScaleBar(width, varargin) adds a scale bar to a figure.
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


function plotScaleBar(width, varargin)

% fixed parameters
b = 40;

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('width', @isscalar);
ip.addParamValue('Handle', gca, @ishandle)
ip.addParamValue('Location', 'southwest', @ischar);
ip.addParamValue('Label', [], @ischar);
ip.addParamValue('FontName', 'Helvetica', @ischar);
ip.addParamValue('FontSize', 10, @isscalar);
ip.parse(width, varargin{:});
label = ip.Results.Label;
fontName = ip.Results.FontName;
fontSize = ip.Results.FontSize;

XLim = get(ip.Results.Handle, 'XLim');
YLim = get(ip.Results.Handle, 'YLim');

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


textProps = {'Color', 'w',...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Center',...
                'FontName', fontName, 'FontSize', fontSize};%, 'FontUnits', 'normalized'};

hold on;
set(gcf, 'InvertHardcopy', 'off');
switch ip.Results.Location
    case 'northeast'
        fill([lx-width lx lx lx-width]-dx, [height height 0 0]+dx,...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(lx-dx-width/2, dx+height, label, textProps{:});
        end
    case 'southeast'
        patch([lx-width lx lx lx-width]-dx, [ly ly ly-height ly-height]-max(dx,textHeight),...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(lx-dx-width/2, ly-max(dx,0.9*textHeight), label, textProps{:});
        end
    case 'southwest'
        fill([0 width width 0]+dx, [ly ly ly-height ly-height]-max(dx,textHeight),...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(dx+width/2, ly-max(dx,0.9*textHeight), label, textProps{:});
        end
    case 'northwest'
        fill([0 width width 0]+dx, [height height 0 0]+dx,...
            [1 1 1]*0.9999, 'EdgeColor', 'none');
        if ~isempty(label)
            text(dx+width/2, dx+height, label, textProps{:});
        end
end
