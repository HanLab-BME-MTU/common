
clear all; clc;

% Data
s = -0.5;
cP = [0 -4 0; 1 3 0; 6 -9 0];

% Test
t = arcLengthToNativeBezierParametrization(cP,s)

% Plot
t1 = linspace(0,1,100)';
ub = max(1,s); lb = min(0,s);
t2 = linspace(lb,ub,100)';
X = renderBezier(cP,t1);
X2 = renderBezier(cP,t2);
P1 = renderBezier(cP,t);
figure(1);
plot(X(:,1),X(:,2),'b','LineWidth',2);
hold on
plot(X2(:,1),X2(:,2),'r');
scatter(P1(1),P1(2));
axis equal
hold off
