function [H,ppSig,ppMu,ljsAll,Htry] = estimateOptimalBandwidth(X,HRange,varargin)
%ESTIMATEOPTIMALBANDWIDTH estimates the optimal mean shift bandwidth for each point
%
%H = estimateOptimalBandwidth(X,HRange)
%H = estimateOptimalBandwidth(X,HRange,...)
%
%   This function implements the bandwidth estimation in [1] to determine,
%   within the specified bandwidth range, what the optimal mean-shift
%   bandwidth is to use for each point in the input dataset. These
%   per-point bandwidths can then be used in variable-bandwidth mean-shift
%   clustering.
%
%   NOTE: This has currently only been developed for and tested with 2D
%   data!!! Also, this only estimates the diagonal elements of the
%   bandwidth matrix (as in the paper).
%
% References:
% [1] D. Comaniciu, "An Algorithm for Data-Driven Bandwidth Selection", IEEE
% T.P.A.M.I, 25 #2 (2003)
%
% Hunter Elliott
% 5/2013

%% ----------------- Input ---------------- %%

ip = inputParser;
ip.addParamValue('nH',20,@isposint);%Number of bandwidths to try within specified range
ip.addParamValue('HTry',[],@(x)(all(x)>0 && all(diff(x))>0));%Optionally input a vector of bandwidths - this overrides the Hrange and nH inputs
ip.addParamValue('w',3,@isposint);%Window size for calculating local jenson-shannon divergence. Higher values will decrease likelihood of spurious estimates, but may require a larger number of test bandwidths to obtain a reliable solution
ip.addParamValue('NumParallel',6,@isposint);%Number of processors to use for parallel computing. Set to 1 to run serially.
ip.addParamValue('flagUseKDTree', 2, @(x) (isempty(x) || (isscalar(x) && ismember(x,0:2))) );%KDtree flag to pass during clustering. See MeanShiftClustering.m for details. Deafult for this function is to use external, since it's faster. Default in sub-function is matlab internal because the external may not exist and may cause errors e.g. under windows7+r2013b
ip.parse(varargin{:});

p = ip.Results;

nH = p.nH;%decreases annoyingness of code.

[n,d] = size(X);

if nH == 1
    warning('ESTIMATEOPTIMALBANDWIDTH:numBandwidths','Only one bandwidth will be tested, JS stability criteria cannot be applied!');
end

if ~isempty(p.flagUseKDTree)
    extraArgs = {'flagUseKDTree',p.flagUseKDTree};
else
    extraArgs = {};
end

%% ------------ Init ------------ %%

if isempty(p.HTry)
    Htry = logspace(log10(HRange(1)),log10(HRange(2)),p.nH);
else
    Htry = p.HTry;
    nH = numel(Htry);
end

%Setup parallel workers if requested
setPoolSize(p.NumParallel)


%% --------- Per-Bandwidth Clustering ----- %%
%Cluster at each bandwidth, and estimate the mean and sigma of the density
%mode underlying each detected cluster at each bandwidth.

nC = zeros(nH,1);
sigEst = cell(nH,1);%Per-cluster estimated sigmas
muEst = cell(nH,1);%Per-cluster estimated mus
ppSig = nan(n,nH,d);%Sigmas estimated for each point at each bandwidth. We assume diagonal covariance matrices.
ppMu = nan(n,nH,d);%Mean for each point at each bandwidth

ptInd = nan(n,nH);

parfor ih = 1:nH
    
    [clInf,ptInd(:,ih),ptTraj] = MeanShiftClustering(X,Htry(ih),'method','standard','minClusterDistance',Htry(ih),extraArgs{:});
    
    nC(ih) = numel(clInf);
    
    sigEst{ih} = nan(nC(ih),d);
    muEst{ih} = nan(nC(ih),d);
    for u = 1:nC(ih)%for each cluster        
        if clInf(u).numPoints > 2*d %make sure we have more points than the # parameters estimated (this should be bandwidth weighted, but this will work for now)

            %Do least-squares estimate of sigma for this cluster, so long as we
            %have enough points

            muEst{ih}(u,:) = clInf(u).ptClusterCenter;%Center of this cluster is assumed to be mean
            Yi = cellfun(@(x)(x(1:end-1,:)),ptTraj(ptInd(:,ih) == u),'Unif',0);%All traj points for this cluster, exlcluding first which is the point itself
            Yi = vertcat(Yi{:});        
            Mi = cellfun(@(x)(diff(x,1,1)),ptTraj(ptInd(:,ih) == u),'Unif',0);
            Mi = vertcat(Mi{:});                
            sig2v = nan(1,d);
            for v = 1:d
                num = sum(Mi(:,v) .* (muEst{ih}(u,v)-Yi(:,v)));
                den = sum(Mi(:,v) .^2);
                sig2v(v) = Htry(ih)^2 * (num/den - 1);                
            end

            %Estimated sigmas which were <= 0 are a good indicator of a bad
            %estimate! Suppress these.
            if all(sig2v>0)
                sigEst{ih}(u,:) = sig2v;
            end
        end        
    end
    
end

%Assign the sigmas to each point in the cluster for per-point sigmas.
%We do this outside the loop to allow parallelization
for ih = 1:nH
    
    for u = 1:nC(ih)
        %NOTE: We should actually not estimate if #pts too low...
        ppSig(ptInd(:,ih) == u,ih,:) = repmat(sigEst{ih}(u,:),[nnz(ptInd(:,ih)==u) 1]);
        ppMu(ptInd(:,ih) == u,ih,:) = repmat(muEst{ih}(u,:),[nnz(ptInd(:,ih)==u) 1]);        
    end
    
end


%% ----------- Get optimal bandwidth for each point ------- %%
%This is done by finding the bandwidths at which the estimated mu and sigma
%were most stable.



ljsAll = nan(n,nH);
H = nan(d,d,n);
He = nan(n,1);
minJS = nan(n,1);
iMinJS = nan(n,1);

for j = 1:n
    currSig = zeros(d,d,nH);%Since we only estimate diagonal cov matrices, put these in full matrix
    for i = 1:d
        currSig(i,i,:) = ppSig(j,:,i);
    end
    if nH > 1
        ljsAll(j,:) = localJensenShannon(squeeze(ppMu(j,:,:)),currSig,p.w);%Calculates difference between estimated distributions at neighboring bandwidths
        [minJS(j),iMinJS(j)] = min(ljsAll(j,:));
        %H(:,:,j) = diag(squeeze(ppSig(j,iMinJS(j),:)));%The optimal bandwidth is that which minizes this local difference, and is therefore most stable.
    else
        iMinJS(j) = 1;%To still allow least-squares estimation without stability criteria in special cases
    end
    H(:,:,j) = currSig(:,:,iMinJS(j));%The optimal bandwidth is that which minizes this local difference, and is therefore most stable.
    He(j) = Htry(iMinJS(j));%And the bandwidth at which this was estimated.
    
end

