% Francois Aguet, last modified: 10/05/2011

function bsplineDemo()

%---------------
% 1-D example
%---------------
x = 1:20;
f = rand(1,20);
d_fx = 0.1;
fx = 1:d_fx:20;

figure;
plot(x, f, 'k.');
hold on;
plot(fx, interpB3SmoothingSpline1D(fx, f), 'b');
plot(fx, interpB3SmoothingSpline1D(fx, f, 0.5), 'r');
legend('Input signal', '\lambda = 0', '\lambda = 0.5');


%---------------
% 2-D example
%---------------

x_t = [1 0 -1 0];
y_t = [0 1 0 -1];
boundary = 'periodic'; % closed curve
% x_t = [0 2 4 6 8];
% y_t = [0 2 4 6 8];
% boundary = 'symmetric'; % open curve

dt = 0.1;
if (strcmp(boundary, 'periodic')==1)
    t = 1:dt:length(x_t)+1-dt; % +1 last -> first node
else
    t = 1:dt:length(x_t)-dt;
end;

% [xi dx d2x] = interpBSpline1D(t, x_t, 3, boundary);
% [yi dy d2y] = interpBSpline1D(t, y_t, 3, boundary);

lambda = 0;
[fx d_fx d2_fx] = interpB3SmoothingSpline1D(t, x_t, lambda, boundary);
[fy d_fy d2_fy] = interpB3SmoothingSpline1D(t, y_t, lambda, boundary);

figure; plot([x_t x_t(1)], [y_t y_t(1)], 'r.-'); axis equal; hold on; plot(fx, fy, 'ko');
axis tight;
legend('Input signal', ['\lambda = ' num2str(lambda)]);
% curveLength =  sum(sqrt(dx.^2 + dy.^2))*dt;
% fprintf('Curve length: %f\n', curveLength);

curvature = abs(d_fx.*d2_fy - d_fy.*d2_fx) ./ (d_fx.^2 + d_fy.^2).^1.5;
figure;
plot(t, curvature, 'co');
hold on;
plot(t, fx, 'k');
plot(t, fy, 'k--');
plot(t, d_fx, 'r');
plot(t, d_fy, 'r--');
plot(t, d2_fx, 'b-');
plot(t, d2_fy, 'b--');
legend('Curvature', 'f_x(t)', 'f_y(t)', 'f_x''(t)', 'f_y''(t)', 'f_x''''(t)', 'f_y''''(t)');

