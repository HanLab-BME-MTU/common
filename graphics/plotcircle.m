%h = plotcircle(center, radius, varargin) plot circles given an array of centers and radii.
%
% Synopsis: h = plotcircle(center, radius)
%           h = plotcircle(center, radius, 'PropertyName', propertyvalue, ...)
%
% Inputs:  
%         center: Nx2 array of centers
%         radius: Nx1 array of radii
%
% Options: 
%    'EdgeColor': Standard plot properties
%    'FaceColor': "
%    'LineStyle': "
%       'Handle': handle to axis object
%
% Output:
%              h: handles to each circle
%
% Francois Aguet, 10/27/2011

function h = plotcircle(center, radius, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('center', @(c) size(c,2)==2);
ip.addRequired('radius', @(r) any(size(r)==1) && all(r>0));
ip.addParamValue('EdgeColor', 'k');
ip.addParamValue('FaceColor', 'none');
ip.addParamValue('LineStyle', '-');
ip.addParamValue('Handle', gca, @ishandle);
ip.parse(center, radius, varargin{:});

center = ip.Results.center;
radius = ip.Results.radius(:);
nc = size(center,1);

if numel(radius)==1 && nc>1
    radius = repmat(radius, [nc 1]);
end

% plot circles
hold on;
axis(ip.Results.Handle, 'equal');
h = arrayfun(@(i) rectangle('Position', [center(i,1)-radius(i) center(i,2)-radius(i) 2*radius(i) 2*radius(i)],...
    'Curvature', [1,1], 'EdgeColor', ip.Results.EdgeColor, 'FaceColor', ip.Results.FaceColor,...
    'LineStyle', ip.Results.LineStyle, 'Parent', ip.Results.Handle), 1:nc,'UniformOutput',false);
h = [h{:}];

plot(0, 0, ip.Results.LineStyle, 'Color', ip.Results.EdgeColor);
