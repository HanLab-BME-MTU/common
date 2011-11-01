%% Test 1
clc; clear all;

n = 1;
m = 50;
x = linspace(0,2 * pi, m);
y = sin(x) + rand(size(x)) * .0;

figure(1);
plot(x,y,'ro');
hold on;

W = ones(m,2);

[P1, t1] = TLSFitBezierWeightedConstraint1([x' y'],W, n);

W([1:10,m-10:m],1) = 0.1;
W([1:10,m-10:m],2) = 0.1;
% W([1:5,m-5:m],1) = 10;
% W([1:5,m-5:m],2) = 10;

[P2, t2] = TLSFitBezierWeightedConstraint1([x' y'],W, n);

C1 = renderBezier(P1, t1);
C2 = renderBezier(P2, t2);

plot(C1(:,1),C1(:,2),'b.');
plot(C2(:,1),C2(:,2),'r.');

axis equal;
hold off;


%% Test 2
clc; clear all;

n = 1;
m = 60;
x = linspace(0,2 * pi, m);
y = 0.5*sin(x) + rand(size(x)) * .5;
data = [x' y'];
% data = data(randperm(size(data,1)),:);
x = data(:,1)';
y = data(:,2)';

figure(1);
plot(x,y,'ro');

tic;
% w = ones(size(data));
w = [1000*ones(size(data,1),1),1000*ones(size(data,1),1)];
% [P3, t3] = TLSFitBezierWeightedConstraint1(data, w, n);
[P3, t3] = TLSFitBezierWeightedConstrainedCP([data data(:,1)], [w w(:,1)], 3, 1/30);
toc;

C3 = renderBezier(P3, sort(t3));
C3_unordered = renderBezier(P3, t3);
% C3 = renderBezier(P3, linspace(0,1,200)');

hold on;

% Plot nodes and control points
plot(C3(:,1),C3(:,2),'g.-');
% plot(P3(:,1),P3(:,2),'bo');

% Plot residuals
for p=1:m
    plot([x(p) C3_unordered(p,1)],[y(p) C3_unordered(p,2)],'b-');
end

axis equal;
hold off;

