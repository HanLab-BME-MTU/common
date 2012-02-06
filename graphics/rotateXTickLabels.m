%ROTATEXTICKLABELS rotates the labels on the x-axis
%
% INPUTS:      ha : axis handle of the plot to be modified
%         'Angle' : rotation in degrees. Default: 45

% Francois Aguet, 22 Feb 2011
% Last modified: 26 July 2011

function rotateXTickLabels(ha, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('ha', @ishandle);
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('AdjustFigure', true, @islogical);
ip.parse(ha, varargin{:});

xa = get(ha, 'XTick');
xla = get(ha, 'XTickLabel');
if ischar(xla) % label is numerical
    xla = arrayfun(@(i) num2str(str2double(xla(i,:))), 1:size(xla,1), 'UniformOutput', false);
end
set(ha, 'XTickLabel', [], 'Units', 'pixels');

fontName = get(ha, 'FontName');
fontSize = get(ha, 'FontSize');

XLim = get(gca, 'XLim');
YLim = get(ha, 'Ylim');
width = diff(XLim);
height = diff(YLim);


% get height of default text bounding box
h = text(0, 0, ' ', 'FontName', fontName, 'FontSize', fontSize);
extent = get(h, 'extent');
pos = get(gca, 'Position');
shift = extent(4)/height*width/pos(3)*pos(4) * sin(ip.Results.Angle*pi/180)/2;
delete(h);


ht = arrayfun(@(k) text(xa(k)-shift, YLim(1)-0.01*height, xla{k},...
    'FontName', fontName, 'FontSize', fontSize,...
    'Units', 'data', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
    'Rotation', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter, 'Parent', ha),...
    1:length(xa), 'UniformOutput', false);

% largest extent (relative to axes units)
extents = cellfun(@(k) get(k, 'extent'), ht, 'UniformOutput', false);
extents = vertcat(extents{:});

pos = get(ha, 'Position');

lmargin = -min(extents(:,1))/width * pos(3); % normalized units in fig. frame

hx = get(ha, 'XLabel');
maxHeight = max(extents(:,4));
if ~strcmpi(get(hx, 'String'), ' ')
    xlheight = get(hx, 'Extent');
    xlheight = xlheight(4);
else
    xlheight = 0;
end
bmargin = (maxHeight+xlheight)/height * pos(4); % data units -> normalized

if ip.Results.AdjustFigure
    hfig = get(ha, 'Parent');
    fpos = get(hfig, 'Position');
    % expand figure window
    
    if lmargin > pos(1)
        fpos(3) = fpos(3) + lmargin-pos(1);
        pos(1) = lmargin;
    end
    fpos(4) = fpos(4) + bmargin-pos(2);
    pos(2) = bmargin;

    set(hfig, 'Position', fpos, 'PaperPositionMode', 'auto');
    set(ha, 'Position', pos);
end

% shift x-label
pos = get(hx, 'Position');
pos(2) = pos(2) - maxHeight;
set(hx, 'Position', pos, 'VerticalAlignment', 'middle');
