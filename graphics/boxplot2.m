%BOXPLOT2 Box plot grouping multiple sets/categories of data, with error bars and SEM.
%
% Inputs:   prm : matrix or cell array of matrices that contain the box properties:
%                 row 1: mean or median
%                 row 2: optional, SEM
%                 row 2/3: 25th percentile, bottom of box
%                 row 3/4: 75th percentile, top of box
%                 row 4/5: optional, bottom whisker
%                 row 5/6: optional, top whisker
% Options:
%  
%     FaceColor : Nx3 matrix of colors, where N is the number of bars or groups
%     EdgeColor : "
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label
%
% Examples:
%
% 1) Simple box plot
% prm = [3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5];
% figure; boxplot2(prm, 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XLabels', arrayfun(@(k) ['S' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 0);
% 
% 
% 2) Multiple groups
% prm = {[3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5],...
%     [3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5]};
% figure; boxplot2(prm, 'BarWidth', 0.8, 'XLabel', 'x label', 'YLabel', 'y label', ...
%     'XLabels', arrayfun(@(k) ['S' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 0);



% Francois Aguet, 22 Feb 2011 (Last modified: 12/08/2011)

function boxplot2(prm, varargin)

if isnumeric(prm)
    prm = {prm};
end

ng = numel(prm); % # of groups
nb = cellfun(@(c) size(c, 2), prm); % # bars in each group

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm');
ip.addParamValue('FaceColor', jet(max(nb)), @(x) size(x,1)==nb || size(x,1)==ng);
ip.addParamValue('EdgeColor', []);
ip.addParamValue('GroupDistance', 0.5, @isscalar);
ip.addParamValue('BorderWidth', [], @isscalar); 
ip.addParamValue('XLabel', [], @ischar);
ip.addParamValue('XLabels', arrayfun(@(k) num2str(k), 1:sum(ng), 'UniformOutput', false), @(x) iscell(x) && (numel(x)==sum(nb)||numel(x)==ng));
ip.addParamValue('YLabel', ' ', @ischar);
ip.addParamValue('YLim', [], @(x) numel(x)==2);
ip.addParamValue('BarWidth', 0.8, @isscalar);
ip.addParamValue('LineWidth', 2, @isscalar);
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('FontName', 'Helvetica', @ischar);
ip.addParamValue('AxisFontSize', 16, @isscalar);
ip.addParamValue('LabelFontSize', 20, @isscalar);
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('X', [], @(x) numel(x)==ng); % cell array of x-coordinates (groups only)
ip.addParamValue('AdjustFigure', true, @islogical);
ip.parse(prm, varargin{:});

faceColor = ip.Results.FaceColor;
nc = size(faceColor,1);
edgeColor = ip.Results.EdgeColor;
ha = ip.Results.Handle;

if isempty(edgeColor)
    edgeColor = zeros(size(faceColor));
end

bw = ip.Results.BarWidth;
dg = ip.Results.GroupDistance; % distance between groups, in bar widths

% x-coords for groups
xa = cell(1,ng);
if isempty(ip.Results.X)
    xa{1} = 1:nb;
    for k = 2:ng
        xa{k} = xa{1} + xa{k-1}(end) + dg;
    end
else
    dx = min(diff(ip.Results.X));
    for k = 1:ng
        w = (nb-1)/2;
        xa{k} = ip.Results.X(k) + (k-1)*dg + (-w:w)*dx/nb;
    end
end

if isempty(ip.Results.BorderWidth)
    border = (bw+dg)/2;
else
    border = ip.Results.BorderWidth;
end


hold on;
for k = 1:ng
    %nb = size(prm{k},2);
    
    xa{k} = (1:nb) + (k-1)*(nb + dg);
    plotSEM = mod(size(prm{k},1),2)==0;
    
    mu = prm{k}(1,:);
    if plotSEM
        p25 = prm{k}(3,:);
        p75 = prm{k}(4,:);
    else
        p25 = prm{k}(2,:);
        p75 = prm{k}(3,:);
    end
    
    % whiskers (plot first to mask bar at '0')
    if plotSEM && size(prm{k},1)>4
        w1 = prm{k}(5,:);
        w2 = prm{k}(6,:);
        plotWhiskers = 1;
    elseif size(prm{k},1)>3
        w1 = prm{k}(4,:);
        w2 = prm{k}(5,:);
        plotWhiskers = 1;
    else
        plotWhiskers = 0;
    end
    
    if plotWhiskers
        he = errorbar(xa{k}, p25, w1-p25, zeros(size(mu)), 'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth);
        setErrorbarStyle(he, 'bottom', ip.Results.ErrorBarWidth);
        he = errorbar(xa{k}, p75, zeros(size(mu)), w2-p75, 'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth);
        setErrorbarStyle(he, 'top', ip.Results.ErrorBarWidth);
    end
    
    % the box
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [p75; p75; p25; p25; p75; p75];
    
    for b = 1:nb
        if nc==nb
            ci = b;
        else
            ci = k;
        end
        patch(xv(:,b), yv(:,b), faceColor(ci,:), 'EdgeColor', edgeColor(ci,:),...
            'LineWidth', ip.Results.LineWidth);
    end
    
    % mean/median line
    line([lb; rb], [mu; mu], 'Color', [0 0 0], 'LineWidth', ip.Results.LineWidth);
    
    % SEM
    if plotSEM
        sigma = prm{k}(2,:);
        he = errorbar(xa{k}, mu, sigma, 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he, 0.15);
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
if ~isempty(ip.Results.YLim);
    set(ha, 'YLim', ip.Results.YLim);
end

% x label
if ~isempty(ip.Results.XLabel)
    xlabel(ip.Results.XLabel, lfont{:});
end

% x labels
if ip.Results.Angle ~= 0
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter,...
        'AdjustFigure', ip.Results.AdjustFigure);
end

% y label
if ~isempty(ip.Results.YLabel)
    ylabel(ip.Results.YLabel, lfont{:});
end
