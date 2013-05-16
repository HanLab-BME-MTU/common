function js = jensenShannonDivergenceSpec(u,s)

r = numel(s);

a = .5*log( abs(mean(s))/ (prod(abs(s))) .^(1/r) );

U = mean(u);

b = (u - U) .^2;

c = sum(s);

js = a + .5*sum(b ./ c);

% [r,d] = size(s);
% 
% 
% a = .5*log( det(mean(s,3))/ (prod(det(s))) .^(1/r) );
% 
% U = mean(u,3);
% 
% b = (u - U) .^2;
% 
% c = sum(s);
% 
% js = a + .5*sum(b ./ c);

