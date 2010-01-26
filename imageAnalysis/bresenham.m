function P = bresenham(p0, p1)
%BRESENHAM computes the integer positions on a line between the
% positions xS and xE
%
% SYNOPSIS P=bresenham(p0,p1)
%
% INPUT p0 : coordinate of line start point
%       p1 : coordinate of line end point
% 
% OUTPUT P : 2xn matrix with the coordinates of all the 
%            integer positions on the line 
%
% Sylvain Berlemont, 2009

dx = p1(2) - p0(2);
dy = p1(1) - p0(1);

incr1 = [1 1];
incr2 = [0 0];

if dx < 0
    dx = -dx;
    incr1(2) = -1;
end

if dy < 0
    dy = -dy;
    incr1(1) = -1;
end

if dy >= dx
    dqr = 2 * dx;
    dqru = dqr - 2 * dy;
    q = dqr - dy;
    l = dy;
    incr2(1) = incr1(1);
else
    dqr = 2 * dy;
    dqru = dqr - 2 * dx;
    q = dqr - dx;
    l = dx;
    incr2(2) = incr1(2);
end

P = zeros(l + 1, 2);
p = p0;

for d = 1:l+1
    P(d, :) = p;
    
    if (q > 0)
        p = p + incr1;
        q = q + dqru;
    else
        p = p + incr2;
        q = q + dqr;
    end
end

end