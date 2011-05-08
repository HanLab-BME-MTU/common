%BARPLOT2 Bar plot grouping multiple sets/categories of data, with error bars.
%
% INPUTS:   prm : cell array of matrices that contain the box properties:
%                 row 1: height
%                 row 2: optional, error bars
%        colors : Nx3 matrix of colors, where N is the number of bars
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label
%
% Example: barplot2({[3 4; 0.5 0.5]}, 'xlabels', {'Class A', 'Class B'});

% Francois Aguet, 18 March 2011
% Last modified: 8 May 2011

function barplot2(prm, varargin)

if isnumeric(prm)
    prm = {prm};
end

ng = length(prm);
nb = cellfun(@(c) size(c, 2), prm);


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm', @iscell);
ip.addOptional('color', []);
ip.addParamValue('EdgeColor', []);
ip.addParamValue('Group', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('GroupLabels', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('GroupDistance', 1, @isscalar);
ip.addParamValue('xlabel', [], @ischar);
ip.addParamValue('xlabels', [], @(x) all(cellfun(@(y) ischar(y), x)));
ip.addParamValue('ylabel', [], @ischar);
ip.addParamValue('BarWidth', 0.8, @isscalar); 
ip.addParamValue('Rotate', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('ErrorBarPosition', 'top',  @(x) strcmpi(x, 'top') | strcmpi(x, 'both'));
ip.parse(prm, varargin{:});
color = ip.Results.color;
ecolor = ip.Results.EdgeColor;


if strcmpi(ip.Results.Group, 'on')
    grouped = true;
else
    grouped = false;
end


if isempty(color)
    if grouped
        color = arrayfun(@(k) hsv2rgb([rand(1,1) ones(1,2)]), 1:ng, 'UniformOutput', false);
    else
        color = arrayfun(@(k) hsv2rgb([rand(nb(k),1) ones(nb(k),2)]), 1:ng, 'UniformOutput', false);
    end
end

if isempty(ecolor)
    if grouped
        ecolor = arrayfun(@(k) [0 0 0], 1:ng, 'UniformOutput', false);
    else
        ecolor = [0 0 0];
    end
end

    
% if size(color,1)==1
%     color = repmat(color, [totalbars 1]);
% end

lfont = {'FontName', 'Helvetica', 'FontSize', 16};

bw = ip.Results.BarWidth;
dg = ip.Results.GroupDistance; % distance between groups, in bar widths

xa = cell(1,ng);

hold on;
for k = 1:ng
    
    xa{k} = (1:nb(k));
    if k > 1
        xa{k} = xa{k} + xa{k-1}(end) + dg;
    end

    height = prm{k}(1,:);   
    
    % bars
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [height; height; zeros(1,nb(k)); zeros(1,nb(k)); height; height];

    if grouped
        arrayfun(@(b) patch(xv(:,b), yv(:,b), color{k}, 'EdgeColor', ecolor{k}, 'LineWidth', 2), 1:nb(k));
    else
        arrayfun(@(b) patch(xv(:,b), yv(:,b), color{b}, 'EdgeColor', ecolor, 'LineWidth', 2), 1:nb(k));
    end
    
    % error bars
    if size(prm{k},1)>1
        %he = errorbar(xa{k}, height, prm{k}(2,:), 'k', 'LineStyle', 'none', 'LineWidth', 2);
        he = errorbar(xa{k}, height, zeros(size(xa{k})), prm{k}(2,:), 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he);
    end

end
hold off;
box off;

groupCenters = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:ng);

xa = [xa{:}];
set(gca, lfont{:}, 'LineWidth', 1.5,...
    'XTick', xa, 'XLim', [0 xa(end)+1]);


width = diff(get(gca, 'XLim'));
height = diff(get(gca, 'YLim'));

% get height of default text bounding box
h = text(0, 0, ' ', lfont{:});
textHeight = get(h, 'extent');
textHeight = textHeight(4);
extent = textHeight/sqrt(2)/2 * width/height;
delete(h);


% x label
if ~isempty(ip.Results.xlabels)
    set(gca, 'XTickLabel', []);
    xlabels = arrayfun(@(k) text(xa(k)-extent,-0.01*height, ip.Results.xlabels{k},...
        'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right',...
        'Rotation', 45, lfont{:}), 1:length(xa));
    
    maxHeight = max(cellfun(@(x) x(4), arrayfun(@(x) get(x, 'extent'), xlabels, 'UniformOutput', false)));
else
    if strcmp(ip.Results.GroupLabels, 'on')
        set(gca, 'XTick', groupCenters, 'XTickLabel', 1:ng);
    else
        set(gca, 'XTickLabel', 1:sum(nb));
    end
    maxHeight = 0;
end

if ~isempty(ip.Results.xlabel)
    hx = xlabel(ip.Results.xlabel, lfont{:});
    
    position = get(hx, 'Position');
    xlabelHeight = get(hx, 'extent');
    xlabelHeight = xlabelHeight(4) - position(2);
    position(2) = position(2) - 0.02*height - maxHeight;
    set(hx, 'Position', position);
else
    xlabelHeight = 0;
end

% set final axis position
if strcmpi(ip.Results.Rotate, 'on')
    total = (height*1.01 + maxHeight + xlabelHeight) * 1.05;
    position = get(gca, 'Position');
    position([2 4]) = [height*0.02+maxHeight+xlabelHeight height]/total;
    set(gca, 'Position', position);
end

% y label
if ~isempty(ip.Results.ylabel)
    ylabel(ip.Results.ylabel, lfont{:});
end




function setErrorbarStyle(he, de)
if nargin<2
    de = 0.2;
end

he = get(he, 'Children');
xd = get(he(2), 'XData');
xd(4:9:end) = xd(1:9:end) - de;
xd(7:9:end) = xd(1:9:end) - de;
xd(5:9:end) = xd(1:9:end) + de;
xd(8:9:end) = xd(1:9:end) + de;
set(he(2), 'XData', xd);
