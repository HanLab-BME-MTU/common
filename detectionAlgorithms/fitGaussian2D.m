function [prmVect, G] = fitGaussian2D(data, prmVect, mode, mask, origin)
%
% function [xp, yp, A, sigma, c] = fitGaussian2D(data, p, mode)
%
% Input: data: 2-D image array
%        p   : [xp, yp, A, sigma, c]
%        mode: specifies which parameters to estimate
%
% Data is assumed to contain a single spot
% Fitting is performed using the Levenberg-Marquardt algorithm.
%
% Francois Aguet, last modified May 2010

if nargin<4
    mask = [];
end
if nargin<5
    if (size(data,1) ~= size(data,2)) || ~mod(size(data,1),2)
        error('Input must be square with odd side length in this mode.');
    end
    w = (size(data,1)-1)/2;
    xa = -w:w;
    ya = xa;
else
    [height width] = size(data);
    xa = (0:width-1) - origin(1);
    ya = (0:height-1) - origin(2);
end


opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

estIdx = false(1,5); % [x y A s c]
estIdx(regexp('xyAsc', ['[' mode ']'])) = true;

[y,x] = ndgrid(ya, xa);
if sum(estIdx)==1 && estIdx(3)==1
    prmVect = A_CF_Gaussian(data, x, y, prmVect);
else
    p = lsqnonlin(@costGaussian, prmVect(estIdx), [], [], opts, data, x, y, prmVect, estIdx, mask);
    prmVect(estIdx) = p;
end

G = gaussian2D(x, y, prmVect);



function [v, J] = costGaussian(p, data, x, y, prmVect, estIdx, mask)
prmVect(estIdx) = p;

[g J] = gaussian2D(x, y, prmVect);
J(:,estIdx==false) = []; % remove unneeded Jacobian components
v = g - data;
v(mask==1) = [];


function [g J] = gaussian2D(x, y, prmVect)

tmp = num2cell(prmVect);
[xp yp A s c] = deal(tmp{:});

r2 = (x-xp).^2+(y-yp).^2;

g_dA = exp(-r2/(2*s^2));
g = A*g_dA;
g_dxp = (x-xp)./s^2.*g;
g_dyp = (y-yp)./s^2.*g;
g_ds = r2/s^3.*g;
g = g + c;

N = numel(x);
g_dc = ones(N,1);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1]) reshape(g_ds, [N 1]) g_dc];


function prmVect = A_CF_Gaussian(data, x, y, prmVect)
r2 = (x-prmVect(1)).^2+(y-prmVect(2)).^2;
g_dA = exp(-r2/(2*prmVect(4)^2));
prmVect(3) = abs(sum(sum((data-prmVect(5)).*g_dA)) / sum(sum(g_dA.^2)));



% Functions that should be incorporated: bi-variate fitting



% Bi-variate
% function [xp yp A sigma c] = bv_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
% x = lsqnonlin(@cost_bv, [xp yp A sigma(1) sigma(2) sigma(3) c], [], [], opts, data);
% xp = x(1);
% yp = x(2);
% A = x(3);
% sigma = [x(4) x(5) x(6)];
% c = x(7);


% function v = cost_bv(p, data)
% xp = p(1);
% yp = p(2);
% A = p(3);
% sigma_x = p(4);
% sigma_y = p(5);
% theta = p(6);
% b = p(7);
% [x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
%
% r2 = ((x-xp)*cos(theta) + (y-yp)*sin(theta)).^2/sigma_x^2/2 + (-(x-xp)*sin(theta) + (y-yp)*cos(theta)).^2/sigma_y^2/2;
% g = A*exp(-r2) + b;
% v = g - data;

% g_dA = exp(-r2/(2*s^2));
% g_db = A*g_dA;
% g_dxp = (x-xp)./s^2.*g_db;
% g_dyp = (y-yp)./s^2.*g_db;
% g_ds = r2/s^3.*g_db;
% g = g_db + b;

%N = numel(data);
%J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1]) reshape(g_ds, [N 1]) reshape(g_db, [N 1])];
% sigma_x = 2;
% sigma_y = 6;
% theta = pi/3;
% [x,y] = ndgrid(0:0.1:nx-1,0:0.1:nx-1);
% %Cx = 1/2; Cxy = 1/8; Cy = 1/4;
% r2 = ((x-xp)*cos(theta) + (y-yp)*sin(theta)).^2/sigma_x + (-(x-xp)*sin(theta) + (y-yp)*cos(theta)).^2/sigma_y;
% g = exp(-r2);