function js = jensenShannonDivergenceSpec(u,s)
%JENSENSHANNONDIVERGENCE jensen shannon divergence for set of normal distributions
% 
% js = jensenShannonDivergenceSpec(mu,sigma2)
%
% This implements the specialized jenson-shannon divergence which gives a
% measure of the net dissimilarity between a of normal distributions.
% This uses formula 17 in [1]
%
% Input:
%   
%       mu - RxD matrix containing the means of R different D-variate
%       distributions
%
%       sigma2 - DxDxR matrix containing the covariance matrices of R
%       different D-variate normal distributions
%
% Output:
%
%   js - scalar specifying divergence between the R distributions. (0 if
%   they are identical)
%
% Hunter Elliott
%

d = size(s,1);
r = size(s,3);

if any(isnan(s(:))) || any(isnan(u(:)))
    %Saves calculation, avoids "matrix is singular" error on NaN sigmas
    js = NaN;
    return
end

detS = arrayfun(@(x)(det(s(:,:,x))),1:r);

a = .5*log( det(mean(s,3)) / (prod(detS)) .^(1/r) );

U = mean(u,1);

b = (u - repmat(U,[r,1]));

c = inv(sum(s,3));

sum1 = 0;
for j = 1:r    
   sum1 = sum1 +  b(j,:) * c * b(j,:)';
end

js = a + .5*sum1;

