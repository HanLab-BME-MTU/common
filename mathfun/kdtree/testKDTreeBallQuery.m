function testKDTreeBallQuery

X = rand(100000,2);
C = [.5 .5];
R = .2;

[idx,d] = KDTreeBallQuery(X,C,R);
assert(all(d{1} <= R));

plot(X(:,1),X(:,2), 'b.');
hold on;
t = 0:pi/100:2*pi;
plot(C(1) + R * cos(t),C(2) + R * sin(t), 'g-');
plot(X(idx{1},1), X(idx{1},2), 'r.');
axis equal;