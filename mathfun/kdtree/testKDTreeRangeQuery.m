function testKDTreeRangeQuery

X = rand(1000000,2);
C = [.5 .5];
L = [.4 .6];

tic;
idx = KDTreeRangeQuery(X,C,L);
toc

plot(X(:,1),X(:,2), 'b.');
hold on;
plot(X(idx{1},1), X(idx{1},2), 'r.');
axis equal;