function testTLSFitBezier(n)

m = 50;
x = linspace(0,2 * pi, m);
y = sin(x) + rand(size(x)) * .05;
plot(x,y,'ro');
tic;
[P1, t1] = TLSFitBezier([x' y'], n);
[P2, t2] = TLSFitBezierFullParameters([x' y'], n);
toc
C1 = renderBezier(P1, t1);
C2 = renderBezier(P2, t2);
hold on;
plot(C1(:,1),C1(:,2),'b-');
plot(C2(:,1),C2(:,2),'r-');
axis equal;

% Compute the dot-product of the first data point and control point
u = [x(1) - P1(1,1), y(1) - P1(1,2)];
v = [P1(2, 1) - P1(1, 1), P1(2, 2) - P1(1, 2)];
dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
fprintf(1, 'dot product at start point (unconstraint) = %E\n', dot);

u = [x(m) - P1(n+1,1), y(m) - P1(n+1,2)];
v = [P1(n, 1) - P1(n+1, 1), P1(n, 2) - P1(n+1, 2)];
dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
fprintf(1, 'dot product at end point (unconstraint) = %E\n', dot);

u = [x(1) - P2(1,1), y(1) - P2(1,2)];
v = [P2(2, 1) - P2(1, 1), P2(2, 2) - P2(1, 2)];
dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
fprintf(1, 'dot product at start point (full param) = %E\n', dot);

u = [x(m) - P2(n+1,1), y(m) - P2(n+1,2)];
v = [P2(n+1, 1) - P2(n, 1), P1(n+1, 2) - P2(n, 2)];
dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
fprintf(1, 'dot product at end point (full param) = %E\n', dot);
