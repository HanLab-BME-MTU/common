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

hx = get(ha, 'XLabel');
xlabelExtent = get(hx, 'Extent');


xa = get(ha, 'XTick');
xla = get(ha, 'XTickLabel');
if ischar(xla) % if singleton labels
    xla = num2cell(xla);
end
set(ha, 'XTickLabel', []);

fontName = get(ha, 'FontName');
fontSize = get(ha, 'FontSize');
pos = get(ha, 'Position');


YLim = get(ha, 'Ylim');
width = diff(get(gca, 'XLim'));
height = diff(YLim);


% get height of default text bounding box
h = text(0, 0, ' ', 'FontName', fontName, 'FontSize', fontSize);
extent = get(h, 'extent');
textHeight = extent(4);
shift = textHeight/sqrt(2)/4 * width/height;
delete(h);


ht = arrayfun(@(k) text(xa(k)-shift, YLim(1)-0.01*height, xla{k},...
    'FontName', fontName, 'FontSize', fontSize,...
    'Units', 'data', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
    'Rotation', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter),...
    1:length(xa), 'UniformOutput', false);

% largest extent
maxHeight = cellfun(@(k) get(k, 'extent'), ht, 'UniformOutput', false);
maxHeight = vertcat(maxHeight{:});
maxHeight = max(maxHeight(:,4));
maxHeightN = maxHeight/height * pos(4); % data units -> normalized

xlh = xlabelExtent(4)/height * pos(4); % x-label height in normalized units

if ip.Results.AdjustFigure
    % space needed: xlh+maxHeightN-pos(2)
    sp = xlh+maxHeightN-pos(2);% - textHeight/height*pos(4);
    hfig = get(ha, 'Parent');
    fpos = get(hfig, 'Position');
    fpos(4) = fpos(4)*(1+sp);
    set(hfig, 'Position', fpos, 'PaperPositionMode', 'auto');
    
    newAxesHeight = pos(4)/(1+sp);
    newTopMargin = (1-pos(2)-pos(4))/(1+sp);
    pos(2) = 1 - newTopMargin - newAxesHeight;
    pos(4) = newAxesHeight;
    set(ha, 'Position', pos);
end

% shift x-label
pos = get(hx, 'Position');
pos(2) = pos(2) - maxHeight;
set(hx, 'Position', pos, 'VerticalAlignment', 'middle');
