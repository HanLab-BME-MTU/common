%BOXPLOT2 Box plot grouping multiple sets/categories of data, with error bars and SEM.
%
% INPUTS:   prm : cell array of matrices that contain the box properties:
%                 row 1: mean or median
%                 row 2: optional, SEM
%                 row 2/3: 25th percentile, bottom of box
%                 row 3/4: 75th percentile, top of box
%                 row 4/5: optional, bottom whisker
%                 row 5/6: optional, top whisker
%        colors : Nx3 matrix of colors, where N is the number of bars
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label

% Francois Aguet, 22 Feb 2011

function boxplot2(prm, colors, xLabels, yLabel)

fsize = 16;

bw = 0.8; % bar width
dx = (1-bw)/2;

dg = 4*dx;

nbars = sum(cellfun(@(x) size(x,2), prm));
if nargin<2 || isempty(colors);
    colors = ones(nbars,3);
    colors(:,1) = rand(1,nbars);
    colors = hsv2rgb(colors);
end
ng = length(prm);

xa = cell(1,ng);

figure('PaperPositionMode', 'auto');
axes('Position', [0.13 0.2 0.775 0.75]);
hold on;
for k = 1:ng
    nb = size(prm{k},2);
    
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
        he = errorbar(xa{k}, p25, p25-w1, zeros(size(mu)), 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he);
        he = errorbar(xa{k}, p75, zeros(size(mu)), w2-p75, 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he);
    end
    
    % the box
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [p75; p75; p25; p25; p75; p75];
    %patch(xv, yv, 'r', 'LineWidth', 2);

    for b = 1:nb
        patch(xv(:,b), yv(:,b), colors(b+(k-1)*nb), 'LineWidth', 2);
    end
    
    
    % mean/median line
    line([lb; rb], [mu; mu], 'Color', [0 0 0], 'LineWidth', 3);
    
    % SEM
    if plotSEM
        sigma = prm{k}(2,:);
        he = errorbar(xa{k}, mu, sigma, 'k', 'LineStyle', 'none', 'LineWidth', 2);
        setErrorbarStyle(he);
    end
end

box off;
xa = [xa{:}];
set(gca, 'FontName', 'Helvetica', 'FontSize', fsize, 'LineWidth', 1.5,...
    'XTick', xa, 'XLim', [0.5-dx/2 xa(end)+0.5+dx/2]);


width = diff(get(gca, 'XLim'));
height = diff(get(gca, 'YLim'));

% get height of default text bounding box
h = text(0, 0, ' ', 'FontName', 'Helvetica', 'FontSize', fsize);
extent = get(h, 'extent');
extent = extent(4)/sqrt(2)/2 * width/height;
delete(h);

% x label
if nargin>2 && length(xLabels)==nbars
    set(gca, 'XTickLabel', []);
    arrayfun(@(k) text(xa(k)-extent,-0.01*height, xLabels{k},...
        'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right',...
        'Rotation', 45, 'FontName', 'Helvetica', 'FontSize', fsize), 1:length(xa));
else
    set(gca, 'XTickLabel', 1:nbars);
end

% y label
if nargin>3 && ~isempty(ylabel)
    ylabel(yLabel, 'FontName', 'Helvetica', 'FontSize', fsize);
end

% print('-depsc2', '-r300', 'test.eps');


function setErrorbarStyle(he)

he = get(he, 'Children');
xd = get(he(2), 'XData');
de = 0.2;
xd(4:9:end) = xd(1:9:end) - de;
xd(7:9:end) = xd(1:9:end) - de;
xd(5:9:end) = xd(1:9:end) + de;
xd(8:9:end) = xd(1:9:end) + de;
set(he(2), 'XData', xd);
