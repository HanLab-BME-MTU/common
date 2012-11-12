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
%               deltaT   - frame rate
%               cluster  - logical - true for cluster the data 
%               nCluster - scalar  - number of clusters
%Output:
%
%       Protrusion.meanTime    - mean time
%       Protrusion.meanConfInt - mean time confidence interval(based on the input alpha)
%       Protrusion.windows     - structure with all measurements per window (see getPersistenceTime for the measurements)
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
nWin                           = length(TS);[Protrusion,Retraction] = getEdgeMotionPersistence(TS,varargin)

Retraction.meanPersTime        = [];
Retraction.meanPersTimeCI      = [];
Retraction.persTimeClusterCI   = [];
Retraction.persTimeClusterIdx  = [];
Retraction.persTimeClusterMean = [];
Retraction.windows(1:nWin)     = struct('limit',[],'PersTime',[],'BlockOut',[],'MaxVeloc',[],'MeanVeloc',[]);

Protrusion.meanPersTime        = [];
Protrusion.meanPersTimeCI      = [];
Protrusion.persTimeClusterCI   = [];
Protrusion.persTimeClusterIdx  = [];
Protrusion.persTimeClusterMean = [];
Protrusion.windows(1:nWin)     = struct('limit',[],'PersTime',[],'BlockOut',[],'MaxVeloc',[],'MeanVeloc',[]);

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


[Protrusion.meanTimeCI,Protrusion.meanTime] = bootStrapMean(outProtTime,alpha,nBoot);
[Retraction.meanTimeCI,Retraction.meanTime] = bootStrapMean(outRetrTime,alpha,nBoot);


if cluster
    
    [Protrusion.timeClusterCI(:,1:nCluster),Protrusion.timeClusterMean(1:nCluster),Protrusion.timeClusterIdx(1:nCluster)] = ...
                                    clusterWindowsVelocity(outProtTime,nBoot,alpha,nCluster);
                                
    [Retraction.timeClusterCI(:,nCluster+1:2*nCluster),Retraction.timeClusterMean(nCluster+1:2*nCluster),Retraction.timeClusterIdx(nCluster+1:2*nCluster)] = ...
                                    clusterWindowsVelocity(outRetrTime,nBoot,alpha,nCluster);
    
end

end%End of main function

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