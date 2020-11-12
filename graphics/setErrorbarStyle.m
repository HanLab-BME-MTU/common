%setErrorbarStyle(he, pos, de) modifies the width of the error bars

% Francois Aguet, Feb 22 2011 (last modif. 01/21/2012)
% Andrew R. Jamieson, Nov. 2016  - function is deprecated, errorbar no
% longer has "Children" allowing editing of style. 
% See: https://www.mathworks.com/matlabcentral/answers/216357-how-can-i-change-the-bar-style-in-a-errorbar-plot

function setErrorbarStyle(he, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('he');
ip.addOptional('de', 0.2, @isscalar);
ip.addParameter('Position', 'both', @(x) any(strcmpi(x, {'both', 'top', 'bottom'})));
ip.parse(he, varargin{:});
de = ip.Results.de;

% he = get(he, 'Children'); % deprecated 
xd = get(he, 'XData');
if strcmpi(ip.Results.Position, 'bottom')
    xd(4:9:end) = xd(1:9:end);
    xd(5:9:end) = xd(1:9:end);
else
    xd(4:9:end) = xd(1:9:end) - de;
    xd(5:9:end) = xd(1:9:end) + de;
end
if strcmpi(ip.Results.Position, 'top')
    xd(7:9:end) = xd(1:9:end);
    xd(8:9:end) = xd(1:9:end);
else
    xd(7:9:end) = xd(1:9:end) - de;
    xd(8:9:end) = xd(1:9:end) + de;
end
set(he, 'XData', xd);