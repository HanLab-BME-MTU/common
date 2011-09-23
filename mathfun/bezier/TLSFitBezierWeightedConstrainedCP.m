function [P t res] = TLSFitBezierWeightedConstrainedCP(data, w, n, varargin)
%
% WORK IN PROGRESS!
%

% TLSFitBezier computes the best total least squares fit of a nth-degree
% Bezier curve to a set of ordered data points. As oppose to the function
% TLSFitBezier, the optimization is performed on the variable t1, ..., tm,
% x0, ..., xn, y0, ... yn.
%
% Required Inputs:
% data           A m x d array representing a set of d-dimensional points
% n              Degree of the Bezier curve.
%
% Optional Inputs:
% MaxFunEvals    Maximum number of fonctional evaluations during lsqnonlin.
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
% Sylvain Berlemont, August 2011

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isnumeric);
ip.addRequired('n', @(n) n > 0);
ip.addParamValue('MaxFunEvals', 1e4, @isscalar);
ip.addParamValue('MaxIter', 1e4, @isscalar);
ip.addParamValue('Display', 'on', @isstring);
ip.addParamValue('TolX', 1e-8, @isscalar);
ip.addParamValue('TolFun', 1e-8, @isscalar);

ip.parse(data, n, varargin{:});
maxFunEvals = ip.Results.MaxFunEvals;
% maxIter = ip.Results.MaxIter;
maxIter = 10;
display = ip.Results.Display;
tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

% Define options of lsqnonlin algorithm
opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', maxFunEvals, ...
    'MaxIter', maxIter, ...
    'Display', display, ...
    'TolX', tolX, ...
    'DerivativeCheck', 'off', ...
    'Tolfun', tolFun);

optsLin = optimset('Jacobian', 'on', ...
    'MaxFunEvals', maxFunEvals, ...
    'MaxIter', maxIter, ...
    'Display', display, ...
    'TolX', tolX, ...
    'DerivativeCheck', 'off', ...
    'LargeScale', 'off', ...
    'Tolfun', tolFun);

% Define initial vector of nodes t0
[m dim] = size(data);
t = linspace(0,1,m)';

% Compute the weights diagonal matrix
W = sparse(diag(w(:)));

% Compute an initial set of control points
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

P(1,:) = data(1,:);
P(n+1,:) = data(end,:);
P(2:n,:) = repmat(P(1,:),n-1,1) + repmat(1/n*(1:n-1)',1,dim) .* repmat((P(n+1,:)-P(1,:)),n-1,1);

resnormOld = -1;
minRad = 5; % Minimal curvature radius

for i=1:maxIter
    % Solve the non-linear optimization on t
    Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
    fun = @(t) r(t, data, W, Cnk, Cn_1k, n, P);
    t = t(2:end-1);
    lb = zeros(size(t));
    ub = ones(size(t));
    [t, ~, res] = lsqnonlin(fun, t, lb, ub, opts);
    t = [0; t; 1];
    
    % Compute Bernstein Matrix
    B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;
    z = zeros(size(B));
    
    % Solve the linear optimization on P
    if n ~= 1
        planeNormal = (P(n+1,:)-P(1,:));
        planePoints = repmat(P(1,:),n-1,1) + repmat(1/n*(1:n-1)',1,dim) .* repmat(planeNormal,n-1,1);
        d = planePoints * planeNormal';
        
        % Degree of freedom reduction constraints
        Aeq_row = [zeros(1,dim); planeNormal; zeros(n-1,dim)];
        Aeq_row = Aeq_row(:)';
        Aeq = Aeq_row;
        for k=1:n-2
            Aeq = [Aeq; circshift(Aeq_row,[0,k])];
        end
        beq = d;
        
        % Pseudo curvature constraints
        alpha = minRad - sqrt(minRad^2-0.25*norm(planeNormal)^2);
        
        ub = inf(n+1,dim);
        ub(2:end-1,:) = planePoints+alpha;
        ub = ub(:);
        lb = -inf(n+1,dim);
        lb(2:end-1,:) = planePoints-alpha;
        lb = lb(:);
        
        % TODO: Not yet 3D
        %     [I J K] = meshgrid([1,-1],[1,-1],[1,-1]);
        %     C = [I(:), J(:), K(:)];
        [I J] = meshgrid([1,-1],[1,-1]);
        C = [I(:), J(:)];
        
        block = zeros(2^dim,(n+1)*dim);
        block(:,2:n+1:dim*(n+1)+1) = C;
        A = block;
        for d=1:n-2
            A = [A; circshift(block,[0,d])];
        end
        
        b = sum(reshape(repmat(planePoints',2^dim,1),dim,2^dim*(n-1))'.*repmat(C,n-1,1),2)+sqrt(2)*alpha;

        C = W*[B z; z B];
        d = W*data(:);
        [P, resnorm, res] = lsqlin(C, d, A, b, Aeq, beq, lb, ub, P(:), optsLin);
    else
        C = W*[B z; z B];
        d = W*data(:);
        [P, resnorm, res] = lsqlin(C, d, [], [], [], [], [], [], P(:), optsLin);
    end
    
    P = reshape(P(1:end), [n+1,dim]);
    
    % ---
    
    % Solve the non-linear optimization on t
    Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
    fun = @(t) r(t, data, W, Cnk, Cn_1k, n, P);
    t = t(2:end-1);
    lb = zeros(size(t));
    ub = ones(size(t));
    [t, ~, res] = lsqnonlin(fun, t, lb, ub, opts);
    t = [0; t; 1];
    
    % ---
    
    if resnorm-resnormOld < tolFun
        disp('Exit criterion satisfied!');
        break;
    else
        resnormOld = resnorm;
    end
end

% Reshape residual
res = sqrt(sum(reshape(res, [m, dim]).^2, 2));

function [F J] = r(t, data, W, Cnk, Cn_1k, n, P)

[m dim] = size(data);

% Append t1 and tm
t = [0; t; 1];

% Compute Bernstein Matrix
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

% Compute the residual
F = data - B * P;
F = W * F(:);

if nargout > 1
    % Compute the Benstein Matrix of order n-1
    Bn_1 = (bsxfun(@power, t, 0:n-1) .* bsxfun(@power, 1 - t, n-1:-1:0)) * Cn_1k;
    
    % Compute derivative of B against t
    z = zeros(m,1);
    Bt = n * ([z Bn_1] - [Bn_1 z]);
    % Bt(1, :) = 0;
    % Bt(end, :) = 0;
    
    % TODO: Only 2D for now! 
    J11 = W(1:m,1:m)*diag(-Bt * P(:,1));
    J21 = W(m+1:end,m+1:end)*diag(-Bt * P(:,2));

    z = zeros(m, n+1);
    J = [J11(:,2:end-1); J21(:,2:end-1)];
end

