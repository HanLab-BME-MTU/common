function TLSFitBezierWeightedConstraint1Test()
clc;
clear all;

n = 1;
m = 50;
x = linspace(0,2 * pi, m);
y = sin(x) + rand(size(x)) * .0;

figure(1);
plot(x,y,'ro');
hold on;

W = ones(m,2);
W([1:10,m-10:m],1) = 0.1;
W([1:10,m-10:m],2) = 0.1;
% W([1:5,m-5:m],1) = 10;
% W([1:5,m-5:m],2) = 10;

[P1, t1] = TLSFitBezierConstraint1([x' y'], n);
[P2, t2] = TLSFitBezierWeightedConstraint1([x' y'],W, n);

C1 = renderBezier(P1, t1);
C2 = renderBezier(P2, t2);

plot(C1(:,1),C1(:,2),'b.');
plot(C2(:,1),C2(:,2),'r.');

axis equal;
hold off;

P2