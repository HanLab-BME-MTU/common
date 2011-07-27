%setErrorbarStyle(he, pos, de) modifies the width of the error bars

% Francois Aguet, Feb 22 2011 (last modif. 07/27/2011)

function setErrorbarStyle(he, pos, de)

if nargin<2
    pos = 'both';
end
if nargin<3
    de = 0.2;
end

he = get(he, 'Children');
xd = get(he(2), 'XData');
if strcmpi(pos, 'bottom')
    xd(4:9:end) = xd(1:9:end);
    xd(5:9:end) = xd(1:9:end);
else
    xd(4:9:end) = xd(1:9:end) - de;
    xd(5:9:end) = xd(1:9:end) + de;
end
if strcmpi(pos, 'top')
    xd(7:9:end) = xd(1:9:end);
    xd(8:9:end) = xd(1:9:end);
else
    xd(7:9:end) = xd(1:9:end) - de;
    xd(8:9:end) = xd(1:9:end) + de;
end
set(he(2), 'XData', xd);