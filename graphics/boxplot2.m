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
%     'XLabels', arrayfun(@(k) ['Group label ' num2str(k)], 1:2, 'UniformOutput', false),...
%     'Angle', 45, 'FaceColor', [1 0.5 0.5; 0.5 1 0.5], 'EdgeColor', [0.8 0 0; 0 0.8 0]);

% Francois Aguet, 22 Feb 2011 (Last modified: 08/15/2012)

function h = boxplot2(prm, varargin)

if isnumeric(prm)
    prm = {prm};
end

ng = size(prm{1},2); % #groups
nb = numel(prm); % #bars in each group

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('prm');
ip.addParamValue('FaceColor', jet(max(nb)), @(x) size(x,1)==1 || size(x,1)==nb || size(x,1)==ng);
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
ip.addParamValue('Interpreter', 'tex', @(x) any(strcmpi(x, {'tex', 'latex', 'none'})));
ip.addParamValue('X', [], @(x) numel(x)==ng); % cell array of x-coordinates (groups only)
ip.addParamValue('AdjustFigure', true, @islogical);
ip.addParamValue('ErrorbarColor', []);
ip.parse(prm, varargin{:});

faceColor = ip.Results.FaceColor;
if size(faceColor,1)==1
    faceColor = repmat(faceColor, [nb 1]);
end
nc = size(faceColor,1);

edgeColor = ip.Results.EdgeColor;
if size(edgeColor,1)==1
    edgeColor = repmat(edgeColor, [nb 1]);
elseif isempty(edgeColor)
    edgeColor = zeros(size(faceColor));
end

errorbarColor = ip.Results.ErrorbarColor;
if isempty(errorbarColor)
    errorbarColor = zeros(size(faceColor));
end

ha = ip.Results.Handle;

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

plotSEM = mod(size(prm{1},1),2)==0;

hold on;
% handles
h = zeros(1,nb);
for k = 1:ng
    
    % concatenate values for group 'k'
    M = cellfun(@(i) i(:,k), prm, 'UniformOutput', false);
    M = [M{:}];
    
    %xa{k} = (1:nb) + (k-1)*(nb + dg);
    
    if plotSEM
        p25 = M(3,:);
        p75 = M(4,:);
    else
        p25 = M(2,:);
        p75 = M(3,:);
    end
    
    % whiskers (plot first to mask bar at '0')
    if plotSEM && size(M,1)>4
        w1 = M(5,:);
        w2 = M(6,:);
        plotWhiskers = 1;
    elseif size(M,1)>3
        w1 = M(4,:);
        w2 = M(5,:);
        plotWhiskers = 1;
    else
        plotWhiskers = 0;
    end
    
    mu = M(1,:);
    
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
        
        if plotWhiskers
            he = errorbar(xa{k}(b), p25(b), w1(b)-p25(b), 0, 'Color', errorbarColor(ci,:), 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'bottom');
            he = errorbar(xa{k}(b), p75(b), 0, w2(b)-p75(b), 'Color', errorbarColor(ci,:), 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
            setErrorbarStyle(he, ip.Results.ErrorBarWidth, 'Position', 'top');
        end
        
        hp = patch(xv(:,b), yv(:,b), faceColor(ci,:), 'EdgeColor', edgeColor(ci,:),...
            'LineWidth', ip.Results.LineWidth);
        if k==1
            h(b) = hp;
        end
        
        % mean/median line
        line([lb(b); rb(b)], [mu(b); mu(b)], 'Color', errorbarColor(ci,:), 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
        
        % replot border
        hp = patch(xv(:,b), yv(:,b), faceColor(ci,:), 'EdgeColor', edgeColor(ci,:),...
            'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
        set(hp, 'FaceColor', 'none');
    end
    
    
    % SEM
    if plotSEM
        sigma = M(2,:);
        he = errorbar(xa{k}, mu, sigma, 'k', 'LineStyle', 'none', 'LineWidth', ip.Results.LineWidth, 'HandleVisibility', 'off');
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

set(ha, 'XTick', la, 'XTickLabel', ip.Results.XLabels, 'XLim', [xa(1)-border xa(end)+border]);
if ~isempty(ip.Results.YLim);
    set(ha, 'YLim', ip.Results.YLim);
end

% x labels
if ip.Results.Angle ~= 0
    rotateXTickLabels(ha, 'Angle', ip.Results.Angle, 'Interpreter', ip.Results.Interpreter,...
        'AdjustFigure', ip.Results.AdjustFigure);
end
