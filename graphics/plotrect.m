function crds = plotrect(rect,opt)
%PLOTRECT plots a rectangle
%
% SYNOPSIS plotrect(rect,opt)
%
% INPUT rect: rectangle [ul(1),ul(2),width,height]
%       opt : plot options (see plot)
%
% OUTPUT crds: 2x4 matrix with corner coordinates 
%
% SEE ALSO plot

x1 = [rect(1),rect(1)+rect(3),rect(1)+rect(3),rect(1),rect(1)];
x2 = [rect(2),rect(2),rect(2)+rect(4),rect(2)+rect(4),rect(2)];

plot(x1,x2,opt);

crds = [x1(1:4);x2(1:4)];
return;