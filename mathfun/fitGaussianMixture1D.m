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

function [mu, sigma, A] = fitGaussianMixture1D(data, n, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isvector);
ip.addRequired('n', @isscalar);
ip.addParamValue('Init', [], @isvector);
ip.addParamValue('Display', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.parse(data, n, varargin{:});
init = ip.Results.Init;

opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-4, ...
    'Tolfun', 1e-4);

[f_ecdf, x_ecdf] = ecdf(data);

if isempty(init)
    mu_init = cumsum(ones(1,n)/(n+1))*2*mean(data);
    sigma_init = std(data)/sqrt(n)*ones(1,n);
    A_init = ones(1,n)/n;
    init = reshape([mu_init; sigma_init; A_init], [1 3*n]);
end
    
[p,resnorm,~,~,~,~,J] = lsqnonlin(@cost, init(1:end-1), [], [], opts, x_ecdf, f_ecdf);

% covariance matrix, error propagation
C = resnorm*full(inv(J'*J));
p_std = sqrt(diag(C)/(numel(data)-length(p) - 1));

% the last Gaussian is constrained; compute variance and append
covA = triu(C(3:3:end,3:3:end));
std_An = sqrt(sum(covA(:))/(numel(data)-length(p) - 1));
p_std = [p_std; std_An]';

mu = [p(1:3:end); p_std(1:3:end)];
sigma = [p(2:3:end); p_std(2:3:end)];
A = [p(3:3:end) 1-sum(p(3:3:end)); p_std(3:3:end)];


if strcmpi(ip.Results.Display, 'on')
    figure; plot(x_ecdf, f_ecdf, 'k', 'LineWidth', 2);
    hold on;
    plot(x_ecdf, model(p, x_ecdf), 'r--', 'LineWidth', 2);
    axis([min(x_ecdf) max(x_ecdf) 0 1]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'YTick', 0:0.1:1);
end


function v = cost(p, x_ecdf, f_ecdf)
v = model(p, x_ecdf) - f_ecdf;


% parameters for each Gaussian: mu, sigma, A
% constraint: sum(A) == 1
function f = model(p, x)
n = (numel(p)+1)/3;
mu = p(1:3:end);
sigma = p(2:3:end);
A = p(3:3:end);
A = [A 1-sum(A)];

f = zeros(size(x));
for k = 1:n
    f = f + A(k)*normcdf(x, mu(k), sigma(k));
end
