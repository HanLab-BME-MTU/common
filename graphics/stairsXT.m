%stairsXT(x, f, varargin) plots a histogram as a single 'step' trace. 

% Francois Aguet, 5/21/2012

function stairsXT(x, f, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('x');
ip.addRequired('f');
ip.addParamValue('EdgeColor', 'k', @(s) ischar(s) || isvector(s));
ip.addParamValue('FaceColor', []);
ip.addParamValue('Bounds', 'closed', @(x) any(strcmpi(x, {'open', 'closed'})));
ip.addParamValue('LineWidth', 1);
ip.parse(x, f, varargin{:});

if numel(x)>2
    dx = x(2)-x(1);
else
    dx = 1;
end
x = [x x(end)+dx];
f = [f f(end)];
if strcmpi(ip.Results.Bounds, 'closed')
    x = [x(1) x x(end)];
    f = [0 f 0];
end

[xb,yb] = stairs(x-0.5, f);
hold on;
if ~isempty(ip.Results.FaceColor)
    patch(xb, yb, ip.Results.FaceColor, 'EdgeColor', 'none')
end
plot(xb, yb, '-', 'Color', ip.Results.EdgeColor, 'LineWidth', ip.Results.LineWidth);
