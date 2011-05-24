function testTLSFitBezier(n)

x = 0:pi/20:2 * pi;
y = sin(x) + rand(size(x)) * .1;
plot(x,y,'ro');
[P, t, res] = TLSFitBezier([x' y'], n);
disp(sqrt(sum(res.^2)));
C = renderBezier(P, t);
hold on;
plot(C(:,1),C(:,2),'b-');
hold off;
