function [Protrusion,Retraction] = getEdgeMotionPersistence(TS,varargin)
% This function bootstrap the mean and confidence intervals for the
% protrusion and retraction persistence time using all time series in input TS
%THIS FUNCTION HAS NO TIME SERIES PRE-PROCESSING
%
%Usage: [Protrusion,Retraction] = getEdgeMotionPersistence(TS,'nBoot',1000,'alpha',0.05)
%
%Input:
%       TS - cell array with one time series in each element
%       Optional:
%               nBoot    - scalar  - number of bootstrap samples
%               alpha    - scalar  - confidence level
%               cluster  - logical - true for cluster the data 
%               nCluster - scalar  - number of clusters
%Output:
%
%       Protrusion.meanTime    - 
%       Protrusion.meanConfInt - 
%       Protrusion.windows     - 
%       
%       Same structure for Retraction
%
%       IF the input "cluster" is TRUE
%           Protrusion.cltDisp  - 
%           Protrusion.cltConfI -
%           Protrusion.index    - cell array with indexes for each cluster 
%
%       Same structure for Retraction    
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
ip.addOptional('cluster',false,@islogical);
ip.addOptional('nCluster',2,@isscalar);
ip.parse(TS,varargin{:});

nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
deltaT   = ip.Results.deltaT;
cluster  = ip.Results.cluster;
nCluster = ip.Results.nCluster;

%**************************************************************************


%%
nWin        = length(TS);
outProtTime = [];
outRetrTime = [];

for iWin = 1:nWin
    
    [Protrusion.windows(iWin),Retraction.windows(iWin)] = getPersistenceTime(TS{iWin},deltaT);
    outProtTime = [outProtTime;Protrusion.windows(iWin).PersTime];
    outRetrTime = [outRetrTime;Retraction.windows(iWin).PersTime];
    
end

% Outcome of this loop is:
%  1 outProtTime -vector with the protrusion time for all time series in TS
%  2 outRetrTime -vector with the retraction time for all time series in TS



%For time series with no Protrusion OR no Retrations - remove NaN
outProtTime(isnan(outProtTime)) = [];
outRetrTime(isnan(outRetrTime)) = [];
%*****************************************************************


[Protrusion.meanConfInt,Protrusion.meanTime] = bootStrapMean(outProtTime,alpha,nBoot);
[Retraction.meanConfInt,Retraction.meanTime] = bootStrapMean(outRetrTime,alpha,nBoot);


if cluster
    
    [Protrusion.cltConfI(:,1:nCluster),Protrusion.cltDisp(1:nCluster),Protrusion.index(1:nCluster)] = ...
                                    clusterWindowsVelocity(outProtTime,nBoot,alpha,nCluster);
                                
    [Retraction.cltConfI(:,nCluster+1:2*nCluster),Retraction.cltDisp(nCluster+1:2*nCluster),Retraction.index(nCluster+1:2*nCluster)] = ...
                                    clusterWindowsVelocity(outRetrTime,nBoot,alpha,nCluster);
else
    
    Retraction.cltConfI = [];
    Retraction.cltDisp  = [];
    Retraction.index    = [];
    
    Protrusion.cltConfI = [];
    Protrusion.cltDisp  = [];
    Protrusion.index    = [];
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