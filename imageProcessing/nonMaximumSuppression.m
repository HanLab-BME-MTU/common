%nms = nonMaximumSuppression(res, th)
%
% Inputs:   res : response
%            th : orientation
%
% Uses the grid conventions of steerableDetector()

% Francois Aguet

function nms = nonMaximumSuppression(res, th)

resXT = padarrayXT(res, [1 1], 'symmetric');

[ny,nx] = size(res);
[x,y] = meshgrid(1:nx,1:ny);
[xi,yi] = meshgrid(0:nx+1, 0:ny+1);

% +1 interp
A1 = interp2(xi, yi, resXT, x+sin(th), y-cos(th));

% -1 interp
A2 = interp2(xi, yi, resXT, x-sin(th), y+cos(th));

nms = res;
nms(res<A1 | res<A2) = 0;