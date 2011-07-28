%BARPLOT2 Bar plot grouping multiple sets/categories of data, with error bars.
%
% INPUTS:   prm : cell array of matrices that contain the box properties:
%                 row 1: height
%                 row 2: optional, error bars
%         color : cell array colors, each containing a Nx3 matrix. Colors cycle through matrix.
%
%       Options : see function content
%
% Examples: 
%
% 1) simple bar plot
% figure; barplot2(rand(1,6), 0.1*rand(1,6), 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XLabels', arrayfun(@(k) ['S' num2str(k)], 1:6, 'UniformOutput', false),...
%     'Angle', 0, 'YLim', [0 1]);
%
% 2) Multiple groups
% figure; barplot2(rand(6,3), 0.1*rand(6,3), 'BarWidth', 1, 'GroupDistance', 1, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XLabels', arrayfun(@(k) ['group ' num2str(k)], 1:6, 'UniformOutput', false),...
%     'Angle', 45, 'YLim', [0 1.2]);


% Note: this function uses patch() since colors can't be controlled with bar()

% Francois Aguet, 18 March 2011 (Last modified: 27 July 2011)

function he = barplot2(prm, varargin)

ng = size(prm,1); % #groups
nb = size(prm,2); % #bars in each group

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm');
ip.addOptional('errorbars', [], @(x) all(size(x)==size(prm)));
ip.addParamValue('Color', jet(nb), @(x) size(x,1)==nb || size(x,1)==ng);
ip.addParamValue('EdgeColor', []);
ip.addParamValue('GroupDistance', 1, @isscalar);
ip.addParamValue('BorderWidth', [], @isscalar); 
ip.addParamValue('XLabel', ' ', @ischar);
ip.addParamValue('XLabels', arrayfun(@(k) num2str(k), 1:sum(ng), 'UniformOutput', false), @(x) iscell(x) && (numel(x)==sum(nb)||numel(x)==ng));
ip.addParamValue('YLabel', ' ', @ischar);
ip.addParamValue('YLim', [], @(x) numel(x)==2);
ip.addParamValue('BarWidth', 0.8, @isscalar); 
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarPosition', 'top',  @(x) strcmpi(x, 'top') | strcmpi(x, 'both'));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('FontName', 'Helvetica', @ischar); % specific
ip.addParamValue('AxisFontSize', 16, @isscalar);
ip.addParamValue('LabelFontSize', 20, @isscalar);
ip.parse(prm, varargin{:});

errorBars = ip.Results.errorbars;
color = ip.Results.Color;
nc = size(color,1);
edgeColor = ip.Results.EdgeColor;
ha = ip.Results.Handle;

if isempty(edgeColor)
    edgeColor = zeros(size(color));
end

bw = ip.Results.BarWidth;
dg = ip.Results.GroupDistance; % distance between groups, in bar widths

if isempty(ip.Results.BorderWidth)
    border = (bw+dg)/2;
else
    border = ip.Results.BorderWidth;
end

xa = cell(1,ng);

hold on;
for k = 1:ng
    
    % x-coords for this group
    xa{k} = 1:nb;
    if k > 1
        xa{k} = xa{k} + xa{k-1}(end) + dg;
    end
    
    height = prm(k,:);
    
    
    % errorbars, if top only
    if ~isempty(errorBars) && strcmpi(ip.Results.ErrorBarPosition, 'top')
        he = errorbar(xa{k}, height, zeros(size(xa{k})), errorBars(k,:),...
            'k', 'LineStyle', 'none', 'LineWidth', 2, 'HandleVisibility', 'off');
        setErrorbarStyle(he, ip.Results.ErrorBarPosition, ip.Results.ErrorBarWidth);
    end
        
    % bars
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [height; height; zeros(1,nb); zeros(1,nb); height; height];

    for b = 1:nb
        if nc==nb
            ci = b;
        else
            ci = k;
        end
        patch(xv(:,b), yv(:,b), color(ci,:), 'EdgeColor', edgeColor(ci,:),...
        'LineWidth', 2);
    end
   
    % errorbars, if two-sided
    if ~isempty(errorBars) && strcmpi(ip.Results.ErrorBarPosition, 'both')
        he = errorbar(xa{k}, height, errorBars(k,:),...
            'k', 'LineStyle', 'none', 'LineWidth', 2, 'HandleVisibility', 'off');
        setErrorbarStyle(he, ip.Results.ErrorBarPosition, ip.Results.ErrorBarWidth);
    end
end
hold off;
box off;

if numel(ip.Results.XLabels)==ng
    la = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:ng);
else
    la = [xa{:}];
end

% position of the bars
xa = [xa{:}];

afont = {'FontName', ip.Results.FontName, 'FontSize', ip.Results.AxisFontSize};
lfont = {'FontName', ip.Results.FontName, 'FontSize', ip.Results.LabelFontSize};

set(ha, afont{:}, 'LineWidth', 1.5,...
    'XTick', la, 'XTickLabel', ip.Results.XLabels, 'XLim', [xa(1)-border xa(end)+border],...
    'TickDir', 'out', 'Layer', 'top');
if ~isempty(ip.Results.YLim)
    set(ha, 'YLim', ip.Results.YLim);
end


% x label
xlabel(ip.Results.XLabel, lfont{:});
ylabel(ip.Results.YLabel, lfont{:});

% x labels
if ip.Results.Angle ~= 0
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle);
end

% re-plot axis on top
% axes('Position', get(ha, 'Position'), 'Box', 'off', 'XTick', [], 'YTick', [],...
%     'HitTest','off', 'Color', 'none', 'LineWidth', 1.5);
