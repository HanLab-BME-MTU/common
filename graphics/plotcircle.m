function plotcircle(ctr,r,opt)
%PLOTCIRCLE plots a circle
%
% SYNOPSIS plotrect(ctr,r,opt)
%
% INPUT ctr: center coordinates [ctr(1),ctr(2)]
%       r :  radius
%       opt : plot options (see plot)
%
% SEE ALSO plot

dt = 2*pi/400;
t = 0:dt:2*pi;

plot(r*cos(t)+ctr(1),r*sin(t)+ctr(2),opt);
axis('equal');
return;
