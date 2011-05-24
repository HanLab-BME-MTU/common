function testTLSFitBezier(n)

x = linspace(0,2 * pi, 50);
y = sin(x) + rand(size(x)) * .5;
plot(x,y,'ro');
tic;
[P, t, res] = TLSFitBezier([x' y'], n);
toc
disp(sqrt(sum(res.^2)));
C = renderBezier(P, t);
hold on;
plot(C(:,1),C(:,2),'b-');
hold off;
