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
% Ouput: hScalebar : handle to the patch graphic object and text graphic
%                    object if applicable
%
% Example: plotScalebar(500, [], 'Label', '500 nm', 'Location', 'SouthEast');
%
%
% Note: there is a bug in '-depsc2' as of Matlab2010b, which misprints patches.
%       When printing, use 'depsc' instead.

% Francois Aguet, March 14 2011 (last modified 07/21/2011)
% Sebastien Besson, July 2011

function hScaleBar=plotScaleBar(width, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('width', @isscalar);
ip.addParamValue('Handle', gca, @ishandle)
ip.addParamValue('Location', 'southwest', @(x) any(strcmpi(x, {'northeast', 'southeast', 'southwest', 'northwest'})));
ip.addParamValue('Label', [], @(x) ischar(x) || isempty(x));
ip.addParamValue('FontName', 'Helvetica', @ischar);
ip.addParamValue('FontSize', [], @isscalar);
ip.addParamValue('Color', [1 1 1], @(x) isvector(x) && numel(x)==3);

ip.parse(width, varargin{:});
label = ip.Results.Label;
fontName = ip.Results.FontName;
fontSize = ip.Results.FontSize;
color = ip.Results.Color;

XLim = get(ip.Results.Handle, 'XLim');
YLim = get(ip.Results.Handle, 'YLim');

lx = diff(XLim);
ly = diff(YLim);

height = width/10; % height of the scale bar
dx = ly/20;

if ~isempty(label)
    if isempty(fontSize)
        fontSize = 1.5*height;        
    end
    % get height of default text bounding box
    h = text(0, 0, label,'FontName', fontName, 'FontSize', fontSize);
    extent = get(h, 'extent');
    textHeight = extent(4) + height;
    textWidth = extent(3);
    delete(h);
else
    textHeight = 0;
    textWidth =0;
end
dx=max(dx,(width+textWidth)/2);

textProps = {'Color', color,...
                'VerticalAlignment', 'Top',...
                'HorizontalAlignment', 'Center',...
                'FontName', fontName, 'FontSize', fontSize};

hold on;
set(gcf, 'InvertHardcopy', 'off');

if ~isempty(strfind(ip.Results.Location,'north')), y0=dx; else y0=ly-height-max(dx,textHeight); end
if ~isempty(strfind(ip.Results.Location,'east')), x0=lx-width-dx; else x0=dx; end

% Create scalebar and optional text
hScaleBar(1) = fill([x0 x0+width x0+width x0],[y0+height y0+height y0 y0],...
    color, 'EdgeColor', 'none');
if ~isempty(label)
    hScaleBar(2) = text(x0+width/2, y0+height, label, textProps{:});
end


