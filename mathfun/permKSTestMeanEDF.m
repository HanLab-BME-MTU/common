function [H, pValue] = permKSTestMeanEDF(d1, d2, varargin)
% permKSTestMeanEDF performs a permutation test between two average EDFs, based on KS-distance
%
% Inputs:
%         d1, d2 : matrices of EDFs from which the mean EDFs to be compared are calculated
%
% Options: 
%          alpha : alpha value, default: 0.05.
%           nrep : # of permutations, default: 1900. This gives a coefficient of
%                  variation <=0.10 for alpha = 0.05. Calculated as
%                  nrep = (1-alpha)/cv^2/alpha. See [1] for details.
%
% Reference:
% [1] Efron, B. and Tibshirani, R., "An introduction to the Bootstrap," Ch. 15, 1993.

% Francois Aguet, 01/14/2013

ip = inputParser;
ip.addRequired('d1', @isnumeric);
ip.addRequired('d2', @isnumeric);
ip.addOptional('alpha', 0.05, @isscalar);
ip.addOptional('nrep', 1900, @isscalar);
ip.parse(d1, d2, varargin{:})
nrep = ip.Results.nrep;

dAll = [d1; d2];

n1 = size(d1,1);
n2 = size(d2,1);
N = n1+n2;

% Calculate the number of permutations. If small, run exact test
w = warning('off', 'MATLAB:nchoosek:LargeCoefficient');
nperms = nchoosek(n1+n2, n1);
warning(w);
if nperms<=nrep % calculate all permutations
    P = false(N, nperms);
    pidx = nchoosek(1:N, n1); % returns row index of class 'sample 1'
    % convert to linear index
    pidx = pidx + repmat(N*(0:nperms-1)', [1 n1]);
    % category (1->sample1, 0->sample2) matrix for all permutations
    P(pidx) = true;
    delta = zeros(nperms,1);
    for i = 1:nperms
        delta(i) = max(abs(mean(dAll(P(:,i),:),1) - mean(dAll(~P(:,i)),1)));
    end
    ns = nperms;
else % compute 'nrep' random permutations
    delta = zeros(nrep,1);
    for i = 1:nrep
        idx = randperm(N); % calculate random permutation of the samples
        delta(i) = max(abs(mean(dAll(idx(1:n1),:),1) - mean(dAll(idx(n1+1:end),:),1)));
    end
    ns = nrep;
end

deltaRef = max(abs(mean(d1,1)-mean(d2,1)));

pValue = sum(delta >= deltaRef)/ns;
H = pValue <= ip.Results.alpha;
