function [meanTime,confI,cltDisp,cltConfI,index,ProtInterval,RetrInterval] = getEdgeMotionPersistence(TS,varargin)
% This function bootstrap the mean and confidence intervals for the
% protrusion and retraction persistence time using all time series in input TS
%
%
%Usage: [meanTime,confI,cltDisp,cltConfI,index] = getEdgeMotionPersistence(TS,'nBoot',1000,'alpha',0.05)
%
%Input:
%       TS - cell array with one time series in each element
%       Optional:
%               nBoot    - scalar  - number of bootstrap samples
%               alpha    - scalar  - confidence level
%               cluster  - logical - true for cluster the data 
%               nCluster - scalar  - number of clusters
%Output:
%       meanTime - 1x2 vector - mean time for the whole data set 
%                  meanTime(1) - for protrusion
%                  meanTime(2) - for retraction
%
%       confI    - 2x2 matrix - confidence interval for the meanTime
%                  confI(1,1) - lower CI for protrusion
%                  confI(2,1) - upper CI for protrusion
%                  confI(1,2) - lower CI for retraction
%                  confI(2,2) - lower CI for retraction
%
%       cltDisp  - vector containing the mean displacement for each cluster
%
%       cltConfI - matrix (2,nCluster) - confidence interval for each
%       cluster
%
%       index    - cell array with indeces for each cluster
%
%See also: getPersistenceTime, findingProtRetrTime
%
%Marco Vilela, 2012

%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('TS',@(x) iscell(x));
ip.addOptional('nBoot',1e3,@isscalar);
ip.addOptional('alpha',.05,@isscalar);
ip.addOptional('deltaT',1,@isscalar);
ip.addOptional('cluster',true,@islogical);
ip.addOptional('nCluster',2,@isscalar);
ip.parse(TS,varargin{:});

nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
deltaT   = ip.Results.deltaT;
cluster  = ip.Results.cluster;
nCluster = ip.Results.nCluster;

%**************************************************************************


%%
nWin = length(TS);
auxP = 1;
auxR = 1;

for iW = 1:nWin
    [ProtTime,RetrTime,ProtInterval{iW},RetrInterval{iW},Up(iW),Dw(iW)] = getPersistenceTime(TS{iW},deltaT);
    outProtTime(auxP:auxP+length(ProtTime) - 1) = ProtTime;
    auxP = auxP + length(ProtTime);
    
    outRetrTime(auxR:auxR+length(RetrTime) - 1) = RetrTime;
    auxR = auxR + length(RetrTime);
    clear ProtTime;clear RetrTime
end

% Outcome of this loop is:
%  1 outProtTime -vector with the protrusion time for all time series in TS
%  2 outRetrTime -vector with the retraction time for all time series in TS



%For time series with no Protrusion OR no Retrations - remove NaN
outProtTime(isnan(outProtTime)) = [];
outRetrTime(isnan(outRetrTime)) = [];
%*****************************************************************


[confI(:,1),meanTime(1)] = bootStrapMean(outProtTime,alpha,nBoot);
[confI(:,2),meanTime(2)] = bootStrapMean(outRetrTime,alpha,nBoot);


if cluster
    [nP,nV] = size(outProtTime);
    if nV > nP
        outProtTime = outProtTime';
    end
    
    [nP,nV] = size(outRetrTime);
    if nV > nP
        outRetrTime = outRetrTime';
    end
    
    [cltConfI(:,1:nCluster),cltDisp(1:nCluster),index(1:nCluster)] = ...
                                    clusterWindowsVelocity(outProtTime,nBoot,alpha,nCluster);
                                
    [cltConfI(:,nCluster+1:2*nCluster),cltDisp(nCluster+1:2*nCluster),index(nCluster+1:2*nCluster)] = ...
                                    clusterWindowsVelocity(outRetrTime,nBoot,alpha,nCluster);
else
    
    cltConfI = [];
    cltDisp  = [];
    index    = [];
    
end

end%End of main function


function [conf,meanS] = bootStrapMean(variable,alpha,nBoot)
%This subFunction bootstrap the mean value of the input "variable" at an
%alpha confidence level with nBoot samples

opt = statset('UseParallel','never');
if matlabpool('size')
    opt = statset('UseParallel','always');
end

[conf,meanSample] = bootci(nBoot,{@nanmean,variable},'alpha',alpha,...
    'type','bca','Options',opt);

meanS = nanmean(meanSample);

end

function [cltConfI,cltDisp,index] = clusterWindowsVelocity(Veloc,nBoot,alpha,nCluster)

    cltDisp  = zeros(1,nCluster);
    cltConfI = zeros(2,nCluster);
    index    = cell(1,nCluster);
    [~,U]    = fcm(Veloc, nCluster);
    maxU     = max(U);
    
    for i =1 : nCluster
        index{i} = find(U(i,:) == maxU);
        [cltConfI(:,i),cltDisp(i)] = bootStrapMean( Veloc(index{i}),alpha,nBoot );
    end
    [cltDisp,idx] = sort(cltDisp); 
    cltConfI      = cltConfI(:,idx); 
end    