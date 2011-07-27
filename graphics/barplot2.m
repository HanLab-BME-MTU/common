%BARPLOT2 Bar plot grouping multiple sets/categories of data, with error bars.
%
% INPUTS:   prm : cell array of matrices that contain the box properties:
%                 row 1: height
%                 row 2: optional, error bars
%         color : cell array colors, each containing a Nx3 matrix. Colors cycle through matrix.
%
%       Options : see function content
%
% Example: 
%
% figure
% barplot2({[3 4 2 3; 0.5 0.5 0.3 0.2], [3 4 2 3; 0.5 0.5 0.3 0.2]},...
%     'xlabels', arrayfun(@(k) ['label ' num2str(k)], 1:8, 'UniformOutput', false),...
%     'xlabel', 'x label', 'ErrorBarPosition', 'top', 'ylabel', 'y label',...
%     'FontSize', 16, 'XLabelFontSize', 20);
% set(gca, 'YLim', [0 5]);

% Francois Aguet, 18 March 2011 (Last modified: 26 July 2011)

function he= barplot2(prm, varargin)

if isnumeric(prm)
    prm = {prm};
end

ng = length(prm); % # of groups
nb = cellfun(@(c) size(c, 2), prm); % # bars in each group

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm', @iscell);
ip.addOptional('color', []);
ip.addParamValue('EdgeColor', []);
ip.addParamValue('GroupDistance', 0.5, @isscalar);
ip.addParamValue('xlabel', [], @ischar);
ip.addParamValue('xlabels', arrayfun(@(k) num2str(k), 1:sum(nb), 'UniformOutput', false), @(x) iscell(x) && (numel(x)==sum(nb)||numel(x)==ng));
ip.addParamValue('ylabel', [], @ischar);
ip.addParamValue('BarWidth', 0.8, @isscalar); 
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarPosition', 'top',  @(x) strcmpi(x, 'top') | strcmpi(x, 'both'));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('FontName', 'Helvetica', @ischar); % specific
ip.addParamValue('FontSize', 16, @isscalar);
ip.addParamValue('XLabelFontSize', 18, @isscalar);
ip.parse(prm, varargin{:});
color = ip.Results.color;
edgeColor = ip.Results.EdgeColor;
ha = ip.Results.Handle;


if isempty(color)
    color = arrayfun(@(k) [0.7 0.9 1], 1:ng, 'UniformOutput', false);
end

if isempty(edgeColor)
    edgeColor = arrayfun(@(k) [0 0.7 1], 1:ng, 'UniformOutput', false);
end


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
    
    % error bars, if top only
    if size(prm{k},1)>1 && strcmpi(ip.Results.ErrorBarPosition, 'top')
        he = errorbar(xa{k}, height, zeros(size(xa{k})), prm{k}(2,:), 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he, ip.Results.ErrorBarWidth, ip.Results.ErrorBarPosition);
    end
        
    % bars
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [height; height; zeros(1,nb(k)); zeros(1,nb(k)); height; height];

    arrayfun(@(b) patch(xv(:,b), yv(:,b), color{k}(mod(b,size(color{k},1))+1,:),...
        'EdgeColor', edgeColor{k}(mod(b,size(edgeColor{k},1))+1,:), 'LineWidth', 2), 1:nb(k));
    
    % error bars, if two-sided
    if size(prm{k},1)>1 && strcmpi(ip.Results.ErrorBarPosition, 'both')
        he = errorbar(xa{k}, height, prm{k}(2,:), 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he, ip.Results.ErrorBarWidth, ip.Results.ErrorBarPosition);
    end
end
hold off;
box off;



if numel(ip.Results.xlabels)==ng
    la = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:ng);
else
    la = [xa{:}];
end

% position of the bars
xa = [xa{:}];

afont = {'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize};
lfont = {'FontName', ip.Results.FontName, 'FontSize', ip.Results.XLabelFontSize};

set(ha, afont{:}, 'LineWidth', 1.5,...
    'XTick', la, 'XTickLabel', ip.Results.xlabels, 'XLim', [0 xa(end)+1],...
    'TickDir', 'out', 'Layer', 'top');


% x label
if ~isempty(ip.Results.xlabel)
    xlabel(ip.Results.xlabel, lfont{:});
end

% x labels
if ip.Results.Angle ~= 0
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle);
end

% y label
if ~isempty(ip.Results.ylabel)
    ylabel(ip.Results.ylabel, lfont{:});
end

% re-plot axis on top
% axes('Position', get(ha, 'Position'), 'Box', 'off', 'XTick', [], 'YTick', [],...
%     'HitTest','off', 'Color', 'none', 'LineWidth', 1.5);



function setErrorbarStyle(he, de, pos)
if nargin<2
    de = 0.2;
end

he = get(he, 'Children');
xd = get(he(2), 'XData');
xd(4:9:end) = xd(1:9:end) - de;
xd(5:9:end) = xd(1:9:end) + de;
if strcmpi(pos, 'top')
    xd(7:9:end) = xd(1:9:end);
    xd(8:9:end) = xd(1:9:end);
else
    xd(7:9:end) = xd(1:9:end) - de;
    xd(8:9:end) = xd(1:9:end) + de;
end
set(he(2), 'XData', xd);
