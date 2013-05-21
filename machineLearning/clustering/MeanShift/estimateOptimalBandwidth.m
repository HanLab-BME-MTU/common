function [H,ppSig,ppMu,ljsAll] = estimateOptimalBandwidth(X,HRange,varargin)
%ESTIMATEOPTIMALBANDWIDTH estimates the optimal mean shift bandwidth for each point
%
%H = estimateOptimalBandwidth(X,HRange)
%
%   This function implements the bandwidth estimation in [1] to determine,
%   within the specified bandwidth range, what the optimal mean-shift
%   bandwidth is to use for each point in the input dataset. These
%   per-point bandwidths can then be used in variable-bandwidth mean-shift
%   clustering.
%
%
%
% References:
% [1] D. Comaniciu, "An Algorithm for Data-Driven Bandwidth Selection", IEEE
% T.P.A.M.I, 25 #2 (2003)
%% ----------------- Input ---------------- %%

ip = inputParser;
ip.addParamValue('nH',20,@isposint);%Number of bandwidths to try within specified range
ip.addParamValue('w',1,@isposint);%Window size for calculating local jenson-shannon divergence.
ip.addParamValue('NumParallel',6,@isposint);

ip.parse(varargin{:});

p = ip.Results;

nH = p.nH;%decreases annoyingness of code.

[n,d] = size(X);


%% ------------ Init ------------ %%

Htry = logspace(log10(HRange(1)),log10(HRange(2)),p.nH);


%Setup parallel workers if requested
if p.NumParallel > 1
    currPoolSize = matlabpool('size');    
    if currPoolSize ~= p.NumParallel       
        if currPoolSize > 0            
            matlabpool('close')
        end        
        matlabpool('open',p.NumParallel);        
    end       
end



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
    
    [clInf,ptInd(:,ih),ptTraj] = MeanShiftClustering(X,Htry(ih),'method','standard','minClusterDistance',Htry(ih));
    
    nC(ih) = numel(clInf);
    
    sigEst{ih} = nan(nC(ih),d);
    muEst{ih} = nan(nC(ih),d);
    for u = 1:nC(ih)%for each cluster
        
        %Do least-squares estimate of sigma for this cluster
        
        muEst{ih}(u,:) = clInf(u).ptClusterCenter;%Center of this cluster is assumed to be mean
        Yi = cellfun(@(x)(x(1:end-1,:)),ptTraj(ptInd(:,ih) == u),'Unif',0);%All traj points for this cluster, exlcluding first which is the point itself
        Yi = vertcat(Yi{:});
        Mi = cellfun(@(x)(diff(x,1,1)),ptTraj(ptInd(:,ih) == u),'Unif',0);
        Mi = vertcat(Mi{:});
        %tu = size(Yi,1);%Number of traj points this cluster            
        
        sig2v = nan(1,d);
        for v = 1:d
            num = sum(Mi(:,v) .* (muEst{ih}(u,v)-Yi(:,v)));
            den = sum(Mi(:,v) .^2);
            sig2v(v) = Htry(ih)^2 * (num/den - 1);                
        end

        sigEst{ih}(u,:) = sqrt(sig2v);
        
        

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

%Since our clustering currently only supports diagonal, isotropic covariane
%matrices, we just use the mean estimated sigma. 
ppSigMean = nanmean(real(ppSig),3);%Temp - Min better? Max? (real is in case bandwidth range was too wide and some estimated sigmas^2 were negative)
ppMuMean = nanmean(real(ppMu),3);

ljsAll = nan(n,nH);
H = nan(n,1);
He = nan(n,1);
minJS = nan(n,1);
iMinJS = nan(n,1);

for j = 1:n
    ljsAll(j,:) = localJensenShannon(ppMuMean(j,:),ppSigMean(j,:),p.w);%Calculates difference between estimated distributions at neighboring bandwidths
    [minJS(j),iMinJS(j)] = min(ljsAll(j,:));
    H(j) = ppSigMean(j,iMinJS(j));%The optimal bandwidth is that which minizes this local difference, and is therefore most stable.
    He(j) = Htry(iMinJS(j));%And the bandwidth at which this was estimated.
    
end

%Clip estimates to input sigma range.
H(H<HRange(1)) = HRange(1);
H(H>HRange(2)) = HRange(2);


