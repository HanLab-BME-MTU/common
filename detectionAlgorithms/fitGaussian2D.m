function [p, G] = fitGaussian2D(data, p, mode)
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
% Francois Aguet, last modified 01/24/2010


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

xp = p(1);
yp = p(2);
A = p(3);
sigma = p(4);
c = p(5);

if (strcmp(mode, 'xyAsc'))
    [xp yp A sigma c] = xyAsc_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'xyAs'))
    [xp yp A sigma] = xyAs_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'xyAc'))
    [xp yp A c] = xyAc_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'xyA'))
    [xp yp A] = xyA_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'xys'))
    [xp yp sigma] = xys_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'xy'))
    [xp yp] = xy_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'As'))
    [A sigma] = As_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts);
elseif (strcmp(mode, 'A'))
    [A] = A_CF_Gaussian(xp, yp, sigma, c, data);
elseif (strcmp(mode, 'bivariate'))
    [xp yp A sigma c] = bv_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts); % sigma: covariance matrix
else
    error('unknown mode');
end;

p = [xp, yp, A, sigma, c];
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
G = A * exp(-r2/(2*sigma^2)) + c;

% figure;
% subplot(1,2,1);
% imagesc(data); colormap(gray(256)); axis image; colorbar;
% subplot(1,2,2);
% imagesc(G); colormap(gray(256)); axis image; colorbar;
% title(['sigma = ' num2str(p(4))]);


function [xp yp A sigma c] = xyAsc_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_xyAsc, [xp yp A sigma c], [], [], opts, data);
xp = x(1);
yp = x(2);
A = x(3);
sigma = x(4);
c = x(5);

function [xp yp A sigma] = xyAs_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_xyAs, [xp yp A sigma], [], [], opts, data, c);
xp = x(1);
yp = x(2);
A = x(3);
sigma = x(4);

function [xp yp A c] = xyAc_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_xyAc, [xp yp A c], [], [], opts, data, sigma);
xp = x(1);
yp = x(2);
A = x(3);
c = x(4);

function [xp yp A] = xyA_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_xyA, [xp yp A], [], [], opts, data, sigma, c);
xp = x(1);
yp = x(2);
A = x(3);

function [xp yp sigma] = xys_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_xys, [xp yp sigma], [], [], opts, data, A, c);
xp = x(1);
yp = x(2);
sigma = x(3);

function [xp yp] = xy_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_xy, [xp yp], [], [], opts, data, sigma, A, c);
xp = x(1);
yp = x(2);

function [A sigma] = As_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_As, [A sigma], [], [], opts, data, xp, yp, c);
A = x(1);
sigma = x(2);

function [A] = A_CF_Gaussian(xp, yp, sigma, c, data)
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*sigma^2));
A = sum(sum((data-c).*g_dA)) / sum(sum(g_dA.^2));

function [xp yp A sigma c] = bv_LvMrq_Gaussian(xp, yp, A, sigma, c, data, opts)
x = lsqnonlin(@cost_bv, [xp yp A sigma(1) sigma(2) sigma(3) c], [], [], opts, data);
xp = x(1);
yp = x(2);
A = x(3);
sigma = [x(4) x(5) x(6)];
c = x(7);



function [v, J] = cost_xyAsc(p, data)
xp = p(1);
yp = p(2);
A = p(3);
s = p(4);
c = p(5);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*s^2));
g = A*g_dA;
g_dxp = (x-xp)./s^2.*g;
g_dyp = (y-yp)./s^2.*g;
g_ds = r2/s^3.*g;
g = g + c;
v = g - data;
N = numel(data);
g_dc = ones(N,1);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1]) reshape(g_ds, [N 1]) g_dc];


function [v, J] = cost_xyAs(p, data, c)
xp = p(1);
yp = p(2);
A = p(3);
s = p(4);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*s^2));
g = A*g_dA;
g_dxp = (x-xp)./s^2.*g;
g_dyp = (y-yp)./s^2.*g;
g_ds = r2/s^3.*g;
g = g + c;
v = g - data;
N = numel(data);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1]) reshape(g_ds, [N 1])];


function [v, J] = cost_xyAc(p, data, s)
xp = p(1);
yp = p(2);
A = p(3);
c = p(4);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*s^2));
g = A*g_dA;
g_dxp = (x-xp)./s^2.*g;
g_dyp = (y-yp)./s^2.*g;
g = g + c;
v = g - data;
N = numel(data);
g_dc = ones(N,1);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1]) g_dc];


function [v, J] = cost_xyA(p, data, s, b)
xp = p(1);
yp = p(2);
A = p(3);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*s^2));
g_db = A*g_dA;
g_dxp = (x-xp)./s^2.*g_db;
g_dyp = (y-yp)./s^2.*g_db;
g = g_db + b;
v = g - data;
N = numel(data);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1])];


function [v, J] = cost_xys(p, data, A, b)
xp = p(1);
yp = p(2);
s = p(3);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_db = A*exp(-r2/(2*s^2));
g_dxp = (x-xp)./s^2.*g_db;
g_dyp = (y-yp)./s^2.*g_db;
g_ds = r2/s^3.*g_db;
g = g_db + b;
v = g - data;
N = numel(data);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_ds, [N 1])];


function [v, J] = cost_xy(p, data, s, A, b)
xp = p(1);
yp = p(2);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*s^2));
g_db = A*g_dA;
g_dxp = (x-xp)./s^2.*g_db;
g_dyp = (y-yp)./s^2.*g_db;
g = g_db + b;
v = g - data;
N = numel(data);
J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1])];


function [v, J] = cost_As(p, data, xp, yp, b)
A = p(1);
s = p(2);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);
r2 = (x-xp).^2+(y-yp).^2;
g_dA = exp(-r2/(2*s^2));
g = A*g_dA;
g_ds = r2/s^3.*g;
g = g + b;
v = g - data;
N = numel(data);
J = [reshape(g_dA, [N 1]) reshape(g_ds, [N 1])];


function v = cost_bv(p, data)
xp = p(1);
yp = p(2);
A = p(3);
sigma_x = p(4);
sigma_y = p(5);
theta = p(6);
b = p(7);
[x,y] = ndgrid(0:size(data,2)-1,0:size(data,1)-1);

r2 = ((x-xp)*cos(theta) + (y-yp)*sin(theta)).^2/sigma_x^2/2 + (-(x-xp)*sin(theta) + (y-yp)*cos(theta)).^2/sigma_y^2/2;
g = A*exp(-r2) + b;

% g_dA = exp(-r2/(2*s^2));
% g_db = A*g_dA;
% g_dxp = (x-xp)./s^2.*g_db;
% g_dyp = (y-yp)./s^2.*g_db;
% g_ds = r2/s^3.*g_db;
% g = g_db + b;
v = g - data;
%N = numel(data);
%J = [reshape(g_dxp, [N 1]) reshape(g_dyp, [N 1]) reshape(g_dA, [N 1]) reshape(g_ds, [N 1]) reshape(g_db, [N 1])];
% sigma_x = 2;
% sigma_y = 6;
% theta = pi/3;
% [x,y] = ndgrid(0:0.1:nx-1,0:0.1:nx-1);
% %Cx = 1/2; Cxy = 1/8; Cy = 1/4;
% r2 = ((x-xp)*cos(theta) + (y-yp)*sin(theta)).^2/sigma_x + (-(x-xp)*sin(theta) + (y-yp)*cos(theta)).^2/sigma_y;
% g = exp(-r2);