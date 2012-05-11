%
% Examples for the binterp() function.
%
% Francois Aguet

%===============================================================================
% Perfect reconstruction for different border conditions
%===============================================================================
f = [4 6 8 10 9 1 2 8];
x = 1:numel(f);

figure;

f_rec = binterp(f, x, 'mirror');
subplot(2,1,1);
plot(x, f-f_rec);

f_rec = binterp(f, x, 'periodic');
subplot(2,1,2);
plot(x, f-f_rec);

%%
%===============================================================================
% 1D example with periodic boundary conditions
%===============================================================================

f = [0 3 2 1 0 3 2 1 0 3 2 1];
x = 1:numel(f);
xi = 1:0.01:numel(f);

[fi, fi_dx, fi_d2x] = binterp(f, xi, bc);

figure;
hold on;
plot(x, f, 'ro');
plot(xi, fi, 'k-');
plot(xi, fi_dx, 'b-');
plot(xi, fi_d2x, 'c-');
legend('f[k]', 'f(x)', 'f''(x)', 'f''''(x)');

%% 
%===============================================================================
% Parametric curves
%===============================================================================
x_t = [1 0 -1 0];
y_t = [0 1 0 -1];

% closed curve
dt = 0.1;
t = 1:dt:length(x_t)+1-dt; % last node == first node

[fx d_fx d2_fx] = binterp(x_t, t, 'periodic');
[fy d_fy d2_fy] = binterp(y_t, t, 'periodic');

figure;
plot([x_t x_t(1)], [y_t y_t(1)], 'r.');
axis equal; hold on;
plot(fx, fy, 'ko');
plotcircle([0 0], 1, 'EdgeColor', 0.6*[1 1 1], 'LineStyle', '--');
axis(1.1*[-1 1 -1 1]);
legend('Samples', 'Cubic B-spline', 'Circle', 'Location', 'NorthEastOutside');

%%
%===============================================================================
% 2D example with periodic boundary conditions
%===============================================================================
f = [0 3 2 1 0 3 2 1 0 3 2 1];
f = f'*f;
[X,Y] = meshgrid(1:0.01:12);

[F F_dx F_dy F_d2x F_d2y] = binterp(f, X, Y, 'periodic');

figure;
subplot(2,2,1)
imagesc(f); colormap(gray(256)); axis image;
title('f');

subplot(2,2,2)
imagesc(F); colormap(gray(256)); axis image;
title('F');

subplot(2,2,3)
imagesc(F_dx); colormap(gray(256)); axis image;
title('d/dx F');

subplot(2,2,4)
imagesc(F_dy); colormap(gray(256)); axis image;
title('d/dy F');

