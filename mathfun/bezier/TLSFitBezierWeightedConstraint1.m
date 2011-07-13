function [P t res] = TLSFitBezierWeightedConstraint1(X, w, n, varargin)
% TLSFitBezier computes the best total least squares fit of a nth-degree
% Bezier curve to a set of ordered data points.
%
% Required Inputs:
% X              A m x d array representing a set of d-dimensional points
% w              A m x d array representing the weight of each point.
%                If the error (standard deviation) of the point is known it
%                can be included by setting w = 1/standard deviation. The
%                full variance/covariance matrix is not supported.
%
% n              Degree of the Bezier curve.
%
% Optional Inputs:
% MaxFunEvals    Aaximum number of fonctional evaluations during lsqnonlin.
% MaxIter        Maximum number of interations during lsqnonlin.
% Display        Verbose mode during lsqnonlin.
% TolX           Tolerance on the solution (i.e. t) during lsqnonlin.
% TolFun         Tolerance on the functional during lsqnonlin.
%
% Outputs:
% P              A n+1 x d array representing the set of d-dimensional
%                control points defining the Bezier curve.
%
% t              A m x 1 array of parameter value for each input points. t
%                belongs to [0, 1].
%
% res            A d*m x 1 array of residue, i.e the components of the distances
%                between a data point and the fitted curve [rx;ry;rz;...]
%
% Sylvain Berlemont, May 2011
% Pascal Berard, July 2011

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('X', @isnumeric);
ip.addRequired('n', @(n) n > 0);
ip.addParamValue('MaxFunEvals', 1e4, @isscalar);
ip.addParamValue('MaxIter', 1e4, @isscalar);
ip.addParamValue('Display', 'off', @isstring);
ip.addParamValue('TolX', 1e-8, @isscalar);
ip.addParamValue('TolFun', 1e-8, @isscalar);

ip.parse(X, n, varargin{:});
maxFunEvals = ip.Results.MaxFunEvals;
maxIter = ip.Results.MaxIter;
display = ip.Results.Display;
tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

% Define options of lsqnonlin algorithm
opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', maxFunEvals, ...
    'MaxIter', maxIter, ...
    'Display', display, ...
    'TolX', tolX, ...
    'Tolfun', tolFun, ...
    'DerivativeCheck','off');

% Define initial vector of nodes t0
[m dim] = size(X);
t = linspace(0,1,m)';
t = t(2:end-1);

% Compute the weights diagonal matrix
WX = w.*X; % Weighted data points
W = zeros(m,m,dim);
for i=1:dim
    W(:,:,i) = diag(w(:,i));
end

% Solve the non-linear optimization on t
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
fun = @(t) r(t, WX, W, Cnk, Cn_1k, n);
[t, ~, res] = lsqnonlin(fun, t, zeros(size(t)), ones(size(t)), opts);

% Compute the control points
t = [0; t; 1];
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

P = zeros(n+1,dim);
for i=1:dim
    WB = W(:,:,i)*B;
    [Q1 R11] = qr(WB,0);
    P(:,i) = R11 \ (Q1' * WX(:,i));
end

% Compute unweighted residuals
res = res./w(:);

end % main function

function [F J] = r(t, WX, W, Cnk, Cn_1k, n)

[m dim] = size(WX);

% Append t1 and tm
t = [0; t; 1];

% Compute Bernstein Matrix
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

% Compute residual = (Id - Bn * pinv(Bn)) * X
WB = zeros(m,n+1,dim);
r = zeros(m,dim);
Q1 = zeros(m,n+1,dim);
R11 = zeros(n+1,n+1,dim);
EE = zeros(1,n+1,dim);
Q2Q2t = zeros(m,m,dim);
for i=1:dim
    WB(:,:,i) = W(:,:,i)*B;
    [Q1(:,:,i) R11(:,:,i) EE(:,:,i)] = qr(WB(:,:,i),0);
    Q2Q2t(:,:,i) = eye(m) - Q1(:,:,i) * Q1(:,:,i)';
    r(:,i) = Q2Q2t(:,:,i) * WX(:,i);
end
F = r(:);

if nargout > 1
    % Compute the Benstein Matrix of order n-1
    Bn_1 = (bsxfun(@power, t, 0:n-1) .* bsxfun(@power, 1 - t, n-1:-1:0)) * Cn_1k;
    
    % Compute derivative of B against t
    z = zeros(m,1);
    Bt = n * ([z Bn_1] - [Bn_1 z]);
    Bt(1, :) = 0;
    Bt(end, :) = 0;
        
    % Compute P    
    WBt = zeros(m,n+1,dim);
    E = zeros(n+1,n+1,dim);
    P = zeros(m,m,dim);
    for i=1:dim
        WBt(:,:,i) = W(:,:,i)*Bt;
        E(sub2ind(size(E), EE(:,:,i), 1:(n+1),repmat(i,1,n+1))) = 1;
        P(:,:,i) = WBt(:,:,i) * E(:,:,i) * (R11(:,:,i) \ Q1(:,:,i)');
    end
        
    % Compute Jacobian Matrix   
    J = zeros(m * dim, m - 2);
    for d = 1:dim
        tmp = -(Q2Q2t(:,:,d) * diag(P(:,:,d) * WX(:,d)) + P(:,:,d)' * diag(r(:,d)));
        J((d-1) * m + 1:d * m,:) = tmp(:, 2:end-1);
    end  
end

end % nested function