%[mu, sigma, A] = fitGaussianMixture1D(data, n, varargin)
%
% Inputs:
%          data : samples of a distribution
%             n : number of Gaussians to fit
%   {'Display'} : 'on' | {'off'} to view results of the fit
%
% Outputs:
%            mu : means of the Gaussians. 2nd row contains propagated error
%         sigma : standard deviations of the Gaussians. 2nd row: propagated error
%             A : amplitudes/relative contributions of the Gaussians. 2nd row: propagated error

% Francois Aguet, 07/19/2011

function [mu, sigma, A, RSS] = fitGaussianMixture1D(data, n, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isvector);
ip.addRequired('n', @isscalar);
ip.addParamValue('Init', [], @isvector);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on','off'})));
ip.addParamValue('Optimizer', 'fmincon', @(x) any(strcmpi(x, {'fmincon','lsqnonlin'})));
ip.addParamValue('ConstrainMeans', false, @islogical);
ip.parse(data, n, varargin{:});
init = ip.Results.Init;

opts = optimset('MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

[f_ecdf, x_ecdf] = ecdf(data);

if isempty(init)
    mu_init = cumsum(ones(1,n)/(n+1))*2*mean(data);
    sigma_init = std(data)/sqrt(n)*ones(1,n);   
    A_init = ones(1,n)/n;
    
    if ip.Results.ConstrainMeans
        init = [mu_init(1) reshape([sigma_init; A_init], [1 2*n])];
        lb = [-Inf zeros(1,2*n)];
        ub = [Inf repmat([Inf 1], [1 n])];
    else
        lb = repmat([-Inf 0 0], [1 n]);
        ub = repmat([Inf Inf 1], [1 n]);
        init = reshape([mu_init; sigma_init; A_init], [1 3*n]);
    end
end


switch ip.Results.Optimizer
    case 'fmincon'
        if ip.Results.ConstrainMeans
            Aeq = [0 repmat([0 1], [1 n])]; % equality constraint
            beq = 1;
            [p,RSS] = fmincon(@(i) sum(costCDFconstr(i, x_ecdf, f_ecdf).^2), init, [], [], Aeq, beq, lb, ub, [], optimset(opts, 'Algorithm', 'interior-point'));
            mu = p(1)*(1:n);
            sigma = p(2:2:end);
            A = p(3:2:end);
        else
            Aeq = repmat([0 0 1], [1 n]); % equality constraint
            beq = 1;
            [p,RSS] = fmincon(@(i) sum(costCDF(i, x_ecdf, f_ecdf).^2), init, [], [], Aeq, beq, lb, ub, [], optimset(opts, 'Algorithm', 'interior-point'));
            mu = p(1:3:end);
            sigma = p(2:3:end);
            A = p(3:3:end);
        end
    case 'lsqnonlin'
        if ip.Results.ConstrainMeans
            [p,RSS,~,~,~,~,J] = lsqnonlin(@costCDFconstr, init, lb, ub, opts, x_ecdf, f_ecdf);
            mu = p(1)*(1:n);
            sigma = p(2:2:end);
            A = p(3:2:end);
        else            
            [p,RSS,~,~,~,~,J] = lsqnonlin(@costCDF, init, lb, ub, opts, x_ecdf, f_ecdf);
            mu = p(1:3:end);
            sigma = p(2:3:end);
            A = p(3:3:end);
        end
end

% % covariance matrix, error propagation
% C = RSS*full(inv(J'*J));
% p_std = sqrt(diag(C)/(numel(data)-length(p) - 1));
%
% % the last Gaussian is constrained; compute variance and append
% covA = triu(C(3:3:end,3:3:end));
% std_An = sqrt(sum(covA(:))/(numel(data)-length(p) - 1));
% p_std = [p_std; std_An]';
%
% mu = [p(1:3:end); p_std(1:3:end)];
% sigma = [p(2:3:end); p_std(2:3:end)];
% A = [p(3:3:end) 1-sum(p(3:3:end)); p_std(3:3:end)];


if strcmpi(ip.Results.Display, 'on')
    figure; plot(x_ecdf, f_ecdf, 'k', 'LineWidth', 2);
    hold on;
    plot(x_ecdf, mixtureModelCDF(x_ecdf, mu, sigma, A), 'r--', 'LineWidth', 2);
    axis([min(x_ecdf) max(x_ecdf) 0 1]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'YTick', 0:0.1:1);
end



function v = costCDF(p, x, f)
mu = p(1:3:end);
sigma = abs(p(2:3:end));
A = abs(p(3:3:end));
v = mixtureModelCDF(x, mu, sigma, A)-f;


function v = costCDFconstr(p, x, f)
n = (numel(p)-1)/2;
mu = p(1)*(1:n);
sigma = abs(p(2:2:end));
A = abs(p(3:2:end));
v = mixtureModelCDF(x, mu, sigma, A)-f;


function f = mixtureModelCDF(x, mu, sigma, A)
f = zeros(size(x));
for i = 1:numel(A)
    f = f + A(i)*normcdf(x, mu(i), sigma(i));
end



