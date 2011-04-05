function testKDTreeRangeQuery

X = rand(100000,2);
C = [.5 .5];
L = [.4 .6];

[idx,d] = KDTreeRangeQuery(X,C,L);

plot(X(:,1),X(:,2), 'b.');
hold on;
t = 0:pi/100:2*pi;
plot(X(idx{1},1), X(idx{1},2), 'r.');
axis equal;