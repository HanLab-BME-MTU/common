function [C,T,res] = snakeBasedBezierFit(data,n,beta,maxNumIter)
% snakeBasedBezierFit computes the optimal n^th Bezier curves that fits
% the data points using a snake-based functional.
%
% data:      a mxd matrix
% n:         an integer > 0
% beta:      a positive scalar
%
% C:         Control points of the optimal Bezier curve. a (n+1)xd matrix
% T:         a mx1 vector
% res:       a mx1 residual vector

[m dim] = size(data);

% problem dimension (number of unknown parameters). Note we do not optimize
% t0 and tm since they need to be t0=0 and tm = 1.
pDim = dim * (n+1) + (m-2);

% Step 1: find an initial guess X = [C, T].
X = zeros(pDim, 1);
% Compute the initial nodes
T = linspace(0,1,m)';
% Compute the initial control points
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
B = (bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0)) * Cnk;
[Q1 R11] = qr(B,0);
C = R11 \ (Q1' * data);
X = [C(:);T(2:end-1)];

dX = ones(pDim,1);
gradF = ones(pDim,1);
numIter = 0;

while norm(dX,2) > 1e-6 && norm(gradF,2) > 1e-6 && numIter < maxNumIter
  
  % Step 3: compute the gradient of F. This is a pDim x 1 vector.
  gradF = computeGradF(X, data, n, beta);
  
  % Step 4: compute the Hessian matrix of F. This is a pDim x pDim matrix.
  hessF = computeHessF(X, data, n, beta);
  
  % Step 5: compute the best direction
  dX = hessF \ gradF;
  
  % TODO: check that gradF . dX < 0 (decreasing direction)
  assert(sum(gradF .* dX) <= eps);
  
  % Step 6: compute the best alpha such that it sufficiently decreases the
  % value of F(X + alpha dX).
  alpha = 1;% computeStepLength(X);
  
  % Step 7: update
  X = X + alpha * dX;
  
  numIter = numIter + 1;
end

% Compute the residual
T = [0; X(dim * (n + 1) + 1:end); 1];
B = (bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0)) * Cnk;
res = sqrt(sum((B * C - data).^2,2));
end

function F = computeF(X, data, n, beta)
[~, dim] = size(data);

% Retrieve the control point coordinates from X
C = repmat(X(1:dim * (n + 1)), n + 1, dim);

% Retrieve the nodes from X and add t0 and tm
T = [0; X(dim * (n + 1) + 1:end); 1];

% Compute the Bernstein matrix
B = bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0);
B = B * diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

% Compute the data fidelity term
dataFidelity = sum(sum((data - B * C).^2, 2));

% Compute the regularization term
regularization = sum(sum(diff(C,1,1),2));

F = dataFidelity + beta * regularization;
end

function gradF = computeGradF(X, data, n, beta)

[m, dim] = size(data);

pDim = dim * (n+1) + (m-2);

% Retrieve the control point coordinates from X
C = reshape(X(1:dim * (n + 1)), n + 1, dim);

% Retrieve the nodes from X and add t0 and tm
T = [0; X(dim * (n + 1) + 1:end); 1];

% Compute the Bernstein matrix
B = bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0);
B = B * diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

%% pre-compute different terms involving factorials and binomials
j = (0:n-1)';
k = 1:n-1;
fact_j = arrayfun(@factorial, j);                        % j!
fact_2n_1 = factorial(2 * n - 1);                        % (2n-1)!
fact_2n_j_2 = arrayfun(@factorial, 2 * n - j - 2);       % (2 * n - j - 2)!
fact_n_j_1 = arrayfun(@factorial, n + j - 1);            % (n + j - 1)!
binom_n_1_j = arrayfun(@(j) nchoosek(n-1,j), j);
binom_n_1_k = arrayfun(@(k) nchoosek(n-1,k), k);
binom_n_1_k_1 = arrayfun(@(k) nchoosek(n-1,k-1), k);

C1 = binom_n_1_j .* fact_j .* fact_2n_j_2 / fact_2n_1;

C2 = fact_n_j_1 ./ fact_j * factorial(n-1) / fact_2n_1;

C3 = bsxfun(@times, binom_n_1_j, bsxfun(@times, 2 * n - ...
  bsxfun(@plus, j, k) - 1, binom_n_1_k_1) - bsxfun(@times, ...
  bsxfun(@plus, j, k), binom_n_1_k)) .* factorial(bsxfun(@plus, j, k) - 1) ...
  .* factorial(2 * n - bsxfun(@plus,j,k) - 2) / fact_2n_1;
%%

% Compute the first term of dF/dxk. dFdC is a (n+1) x dim matrix where
% dFdC(k,l) = 1st term of the derivative of F with respect to the lth
% coordinate of the kth control point.
dFdC1 = 2 * B' * (B * C - data);

% Compute the second term of dF/dxk with k = 0. dFdC2_k0 is a 1 x dim
% vector where dFdC2_k0(l) = 2nd term of the derivative of F with respect
% to the lth coordinate of the first control point (k=0).
dFdC2_k0 = -2 * beta * n^2 * sum(bsxfun(@times, diff(C,1,1), C1), 1);

% Compute the second term of dF/dxk with k = n. dFdC2_kn is a 1 x dim
% vector where dFdC2_kn(l) = 2nd term of the derivative of F with respect
% to the lth coordinate of the last control point (k=n).
dFdC2_kn = 2 * beta * n^2 * sum(bsxfun(@times, diff(C,1,1), C2), 1);

% Compute the second term of dF/dxk with 0 < k < n. dFdC2_k is (n-1) x dim
% matrix where dFdC2_k(k,l) = 2nd term of the derivative of F with respect
% to the lth coordinate of the kth point.
dFdC2_k = 2 * beta * n^2 * permute(sum(bsxfun(@times, permute(C3, [1 3 2]), ...
  diff(C,1,1)), 1), [3 2 1]);

gradF = zeros(pDim, 1);

gradF(1:dim * (n + 1)) = dFdC1(:) + reshape([dFdC2_k0; dFdC2_k; dFdC2_kn], ...
  dim * (n + 1), 1);

% Compute dF/dtk
gradF(dim * (n + 1) + 1:end) = 2 * n * sum((B(2:end-1,:) * C - ...
  data(2:end-1,:)) .* (B(2:end-1,1:end-1) * diff(C,1,1)),2);

end

function hessF = computeHessF(X, data, n, beta)

[m, dim] = size(data);

pDim = dim * (n+1) + (m-2);

% Retrieve the control point coordinates from X
C = reshape(X(1:dim * (n + 1)), n + 1, dim);

% Retrieve the nodes from X and add t0 and tm
T = [0; X(dim * (n + 1) + 1:end); 1];

% Compute the Bernstein matrix
B = bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0);
B = B * diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

%% pre-compute different terms involving factorials and binomials
k = 1:n-1;
l = k';
fact_k = arrayfun(@factorial, k);                        % k!
fact_k_1 = arrayfun(@factorial, k - 1);                  % (k-1)!
fact_n_1 = factorial(n - 1);                             % (n-1)!
fact_2n_1 = factorial(2 * n - 1);                        % (2n-1)!
fact_2n_k_1 = arrayfun(@factorial, 2 * n - k - 1);       % (2 * n - k - 1)!
fact_2n_k_2 = arrayfun(@factorial, 2 * n - k - 2);       % (2 * n - k - 2)!
fact_n_k_2 = arrayfun(@factorial, n + k - 2);            % (n + k - 2)!
binom_n_k = arrayfun(@(k) nchoosek(n,k), k);
binom_n_1_k = arrayfun(@(k) nchoosek(n-1,k), k);
binom_n_1_k_1 = arrayfun(@(k) nchoosek(n-1,k-1), k);

C1 = (binom_n_1_k_1 .* fact_k_1 .* fact_2n_k_1 - binom_n_1_k .* fact_k .* ...
  fact_2n_k_2) / fact_2n_1;
C2 = fact_n_1 * fact_n_k_2 ./ (fact_k_1 * fact_2n_1) .* (1 - (n+k-1) ./ k);
C3 = (bsxfun(@times, bsxfun(@plus, k, l - 1), binom_n_1_k) .* bsxfun(@plus, ...
  (1 - 2 * n) * binom_n_1_k_1', bsxfun(@times, bsxfun(@plus,k,l), binom_n_k')) ...
  + bsxfun(@times, bsxfun(@plus, k, l + 1 - 2 * n), binom_n_1_k_1) .* ...
  bsxfun(@plus, - 2 * n * binom_n_1_k_1' - binom_n_1_k', bsxfun(@times, ...
  bsxfun(@plus, k, l), binom_n_k'))) .* factorial(bsxfun(@plus, k, l - 2)) ...
  .* factorial(bsxfun(@plus, 2 * n - k, -l - 2)) / fact_2n_1;
%%

hessF = zeros(pDim);

% Compute d2Fdxkdxl

% k = 0, l = 0
hessF(1,1) = 2 * sum(B(:,1).^2) + 2 * beta * n^2 / (2 * n - 1);

% k = 0, l = n and symetric case
hessF(1,n+1) = 2 * sum(B(:,1) .* B(:,end)) - 2 * beta * n^2 * fact_n_1^2 / fact_2n_1;
hessF(n+1,1) = hessF(1,n+1);

% k = n, l = n
hessF(n+1,n+1) = 2 * sum(B(:,end).^2) + 2 * beta * n^2 / (2 * n - 1);

% k = 0, 0 < l < n and symetric case
hessF(1,2:n) = 2 * sum(bsxfun(@times, B(:,1), B(:,2:n)), 1) - 2 * beta * n^2 * C1;
hessF(2:n,1) = hessF(1,2:n);

% k = n, 0 < l < n and symetric case
hessF(n+1, 2:n) = 2 * sum(bsxfun(@times, B(:,end), B(:,2:n)), 1) + 2 * beta * n^2 * C2;
hessF(2:n, n+1) = hessF(n+1, 2:n);

% 0 < k < n, 0 < l < n
hessF(2:n, 2:n) = 2 * B(:,2:end-1)' * B(:,2:end-1) + 2 * beta * n^2 * C3;

blks = cell(dim,1);
[blks{:}] = deal(hessF(1:n+1,1:n+1));
hessF(1:dim * (n+1), 1:dim * (n+1)) = blkdiag(blks{:});

% Compute d2Fdtdx
for iDim = 1:dim
  
  offset = (iDim - 1) * (n + 1) + 1;
  
  % d2Fdtdx is a m x (n+1) matrix
  d2Fdtdx = 2 * B(2:end-1,:) .* (sum(bsxfun(@times, bsxfun(@rdivide, ...
    bsxfun(@minus, bsxfun(@plus, permute(0:n, [3, 1, 2]), 0:n), ...
    2 * n * T(2:end-1)), T(2:end-1) .* (1-T(2:end-1))), permute(...
    bsxfun(@times, B(2:end-1,:), C(:,iDim)'), [1, 3, 2])), 3) - ...
    bsxfun(@times, bsxfun(@rdivide, bsxfun(@minus, 0:n, n * T(2:end-1)), ...
    T(2:end-1) .* (1-T(2:end-1))), data(2:end-1,iDim)));
  
  hessF(end-m+3:end, offset:offset+n) = d2Fdtdx;
  hessF(offset:offset+n, end-m+3:end) = d2Fdtdx';
end

% Compute d2Fdt2
hessF(dim * (n + 1) + 1:end, dim * (n + 1) + 1:end) = diag(...
  sum((B(2:end-1,1:end-1) * diff(C,1,1)).^2, 2) + n * (n-1) * ...
  sum((B(2:end-1,:) * C - data(2:end-1,:)) .* (B(2:end-1,1:end-2) * ...
  diff(C,2,1)),2));
end

function alpha = computeStepLength(X)
end