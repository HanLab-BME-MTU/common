function [P t res] = TLSFitBezier(X, n, varargin)
% TLSFitBezier computes the best total least squares fit of a nth-degree
% Bezier curve to a set of ordered data points.
%
% Required Inputs:
% X              A m x d array representing a set of d-dimensional points
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
% res            A m x 1 array of residue, i.e the orthogonal distance
%                between a data point and the fitted curve.
%
% Sylvain Berlemont, May 2011

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
opts = optimset('Jacobian', 'off', ...
  'MaxFunEvals', maxFunEvals, ...
  'MaxIter', maxIter, ...
  'Display', display, ...
  'TolX', tolX, ...
  'Tolfun', tolFun);

% Define initial vector of nodes t0
m = size(X,1);
t = linspace(0,1,m)';

% Solve the non-linear optimization on t
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
fun = @(t) r(t, X, Cnk, n);
[t, ~, res] = lsqnonlin(fun, t, zeros(size(t)), ones(size(t)), opts);

% Compute the control points
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;
x = B \ X(:,1);
y = B \ X(:,2);
P = [x y];

function F = r(t, X, Cnk, n)
    
% Compute Bernstein Matrix
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

% Compute the QR decomposition of B
% Q1 is a m x (n+1) matrix
%Q1 = qr(B,0);

% Q2 is a m x (m - n + 1) matrix
% Q2 * Q2' is a m x m matrix
%Q2Q2t = eye(numel(t)) - Q1 * Q1';

PP = eye(numel(t)) - B * pinv(B);

F = sqrt(sum((PP * X).^2,2));

