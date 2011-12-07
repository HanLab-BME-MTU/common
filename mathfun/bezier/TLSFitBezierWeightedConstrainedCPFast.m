function [P t res] = TLSFitBezierWeightedConstrainedCPFast(data, w, n, maxCurvature, varargin)
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

% Define options of lsqnonlin algorithm
ip.parse(data, n, varargin{:});
% % % maxFunEvals = ip.Results.MaxFunEvals;
% % % maxIter = ip.Results.MaxIter;
maxIter = 100;
% % % % display = ip.Results.Display;
% % % display = 'off';
% % % tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

[m dim] = size(data);
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

% Fit a line
% TODO: Replace with a weighted TLS fit
[x0,a,~,~] = ls3dline(data); 
[minPoint,maxPoint,~,~,~] = projectionPointsLineMinMax(data,x0',a');

% Compute the planar constraints
planeNormal = (maxPoint-minPoint);
planePoints = repmat(minPoint,n+1,1) + repmat(1/n*(0:n)',1,dim) .* repmat(planeNormal,n+1,1);
d = planePoints * planeNormal';

Aeq_row = [planeNormal; zeros(n,dim)];
Aeq_row = Aeq_row(:)';
Aeq = zeros(n+1,(n+1)*dim);
Aeq(1,:) = Aeq_row;
for k=1:n
    Aeq(k+1,:) = circshift(Aeq_row,[0,k]);
end
beq = d;

% Define initial node vector t
t = linspace(0,1,m)';

% Initial control points P
P = planePoints;

% Compute the weights diagonal matrix
W = sparse(diag(w(:)));

% Minimal curvature radius
minRad = 1/maxCurvature; 

% Value of the objective function
resnormOld = -1;

for i=1:maxIter
    % Solve the non-linear optimization on t
    Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
    [~,t] = levmar('lsqnonlinFun','lsqnonlinJac',t,W*data(:),100,[],'unc',m,dim,W,Cnk,Cn_1k,n,P);
    
    % Compute Bernstein Matrix
    B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;
    z = zeros(size(B));
    
    % Solve the linear optimization on P
    if n ~= 1
        % Deviation constraints
        len = norm(planeNormal);        
        if len < 2*minRad
            alpha = minRad - sqrt(minRad^2-0.25*len^2);
        else
            alpha =  len/2;
        end
        
        ub = planePoints+alpha;
        ub = ub(:);
        lb = planePoints-alpha;
        lb = lb(:);
        
        % TODO: Generalize to n dimensions
        if dim == 2
            [I J] = meshgrid([1,-1],[1,-1]);
            C = [I(:), J(:)];
        elseif dim == 3
            [I J K] = meshgrid([1,-1],[1,-1],[1,-1]);
            C = [I(:), J(:), K(:)];
        end
        
        block = zeros(2^dim,(n+1)*dim);
        block(:,1:n+1:dim*(n+1)) = C;
        A = zeros((n+1)*2^dim,(n+1)*dim);
        A(1:2^dim,:) = block;
        for d=1:n
            A(d*2^dim+1:(d+1)*2^dim,:) = circshift(block,[0,d]);
        end
        
        b = sum(reshape(repmat(planePoints',2^dim,1),dim,2^dim*(n+1))'.*repmat(C,n+1,1),2)+sqrt(2)*alpha;
        
        % TODO: Generalize to n dimensions
        if dim == 2
            C = W*[B z; z B];
        elseif dim == 3
            C = W*[B z z; z B z; z z B];
        end
        
        d = W*data(:);        
        [~,P] = levmar('lsqlinFun','lsqlinJac',P(:),d,100,[],'bleic',lb,ub,Aeq,beq,-A,-b,C);
        resnorm = (C*P-d)*(C*P-d)';

    else
        % TODO: Generalize to n dimensions
        if dim == 2
            C = W*[B z; z B];
        elseif dim == 3
            C = W*[B z z; z B z; z z B];
        end
        
        d = W*data(:);       
        [~,P] = levmar('lsqlinFun','lsqlinJac',P(:),d,100,[],'bleic',lb,ub,Aeq,beq,-A,-b,C);
        resnorm = (C*P-d)*(C*P-d)';
    end
    
    P = reshape(P(1:end), [n+1,dim]);

    if resnorm-resnormOld < tolFun
        break;
    else
        resnormOld = resnorm;
    end
end

% Solve the non-linear optimization on t
Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
[~,t] = levmar('lsqnonlinFun','lsqnonlinJac',t,W*data(:),100,[],'unc',m,dim,W,Cnk,Cn_1k,n,P);
res = lsqnonlinFun(t,m,dim,W,Cnk,Cn_1k,n,P)-W*data(:);   
    
% Truncate Bezier curve
[P,t] = truncateBezier(P,min(t),max(t),t);

% Compute unweighted residuals
res = res./w(:);

% Reshape residual
res = reshape(res,[m, dim]);
