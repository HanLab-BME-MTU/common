function testKDTreeBallQuery

X = rand(5000,2);
C = [.5 .5];
R = .2;

[idx,d] = KDTreeBallQuery(X,C,R);

plot(X(:,1),X(:,2), 'b.');
hold on;
ezplot(@(x,y) sqrt((x-C(1)).^2 + (y-C(2)).^2 - R^2), [0, 1, 0, 1]);
plot(X(idx{1},1), X(idx{1},2), 'r.');
axis equal;