%BOXPLOT2 Box plot grouping multiple sets/categories of data, with error bars and SEM.
%
% INPUTS:   prm : matrix or cell array of matrices that contain the box properties:
%                 row 1: mean or median
%                 row 2: optional, SEM
%                 row 2/3: 25th percentile, bottom of box
%                 row 3/4: 75th percentile, top of box
%                 row 4/5: optional, bottom whisker
%                 row 5/6: optional, top whisker
%         color : Nx3 matrix of colors, where N is the number of bars
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label
%
% Example: boxplot2({[3 4; 0.2 0.2; 2 3; 4 5; 0.5 0.5; 0.5 0.5]});

% Francois Aguet, 22 Feb 2011 (Last modified: 27 July 2011)

function boxplot2(prm, varargin)

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
ip.addParamValue('BorderWidth', 0.8, @isscalar); 
ip.addParamValue('Angle', 45, @(x) isscalar(x) && (0<=x && x<=90));
ip.addParamValue('ErrorBarWidth', 0.2, @(x) 0<x && x<=1);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('FontName', 'Helvetica', @ischar); % specific
ip.addParamValue('FontSize', 16, @isscalar);
ip.addParamValue('XLabelFontSize', 18, @isscalar);
ip.addParamValue('YLim', [], @(x) numel(x)==2);
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
        he = errorbar(xa{k}, p25, w1, zeros(size(mu)), 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he, 'bottom', ip.Results.ErrorBarWidth);
        he = errorbar(xa{k}, p75, zeros(size(mu)), w2, 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he, 'top', ip.Results.ErrorBarWidth);
    end
    
    % the box
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [p75; p75; p25; p25; p75; p75];
    
    arrayfun(@(b) patch(xv(:,b), yv(:,b), color{k}(mod(b,size(color{k},1))+1,:),...
        'EdgeColor', edgeColor{k}(mod(b,size(edgeColor{k},1))+1,:), 'LineWidth', 2), 1:nb(k));
    
    % mean/median line
    line([lb; rb], [mu; mu], 'Color', [0 0 0], 'LineWidth', 3);
    
    % SEM
    if plotSEM
        sigma = prm{k}(2,:);
        he = errorbar(xa{k}, mu, sigma, 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he, 0.15);
    end
end
hold off;
box off;

if numel(ip.Results.xlabels)==ng
    la = arrayfun(@(k) (xa{k}(1) + xa{k}(end))/2, 1:ng);
else
    la = [xa{:}];
end

xa = [xa{:}];

afont = {'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize};
lfont = {'FontName', ip.Results.FontName, 'FontSize', ip.Results.XLabelFontSize};

set(ha, afont{:}, 'LineWidth', 1.5,...
    'XTick', la, 'XTickLabel', ip.Results.xlabels, 'XLim', [1-ip.Results.BorderWidth xa(end)+ip.Results.BorderWidth],...
    'TickDir', 'out', 'Layer', 'top');
if ~isempty(ip.Results.YLim);
    set(ha, 'YLim', ip.Results.YLim);
end

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
