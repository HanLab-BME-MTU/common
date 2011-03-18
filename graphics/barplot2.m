%BOXPLOT2 Box plot grouping multiple sets/categories of data, with error bars and SEM.
%
% INPUTS:   prm : cell array of matrices that contain the box properties:
%                 row 1: height
%                 row 2: optional, error bars
%        colors : Nx3 matrix of colors, where N is the number of bars
%       xLabels : cell array of strings, labels for each bar
%        yLabel : string, y-axis label
%
% Example: barplot2({[3 4; 0.5 0.5]}, [], {'Class A', 'Class B'});

% Francois Aguet, 18 March 2011

function barplot2(prm, colors, xLabels, yLabel)

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

    height = prm{k}(1,:);   
    
    % bars
    lb = xa{k} - bw/2;
    rb = xa{k} + bw/2;
    xv = [lb; rb; rb; lb; lb; rb];
    yv = [height; height; zeros(1,nb); zeros(1,nb); height; height];

    for b = 1:nb
        patch(xv(:,b), yv(:,b), colors(b+(k-1)*nb,:), 'LineWidth', 2);
    end
    
    % error bars
    if size(prm{k},1)>1
        he = errorbar(xa{k}, height, prm{k}(2,:), 'k', 'LineStyle', 'none', 'LineWidth', 2);
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

YLim = get(gca, 'YLim');
set(gca, 'Ylim', [0, YLim(2)+1]);


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
