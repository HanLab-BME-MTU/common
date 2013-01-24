function [newTsMap,displayTS,smoothedTsMap,smoothedSNRMap] = getSNRmap(tsMap,varargin)
% This function estimates the local SNR in a signal
%
%Usage:
%       [newTsMap,Up,Dw,snrMap,globalSNR,smoothedTsMap,smoothedSNRMap] = getSNRmap(tsMap,varargin)
%
%Input:
%       tsMap    - matrix with observations_X_variables (Ex:protrusion map)
%       alpha    - alpha level for testing the IMF's (see testImf)
%       outLevel - Std for the outliear detection    (see detectOutliers)
%       trend    - type of trend to be removed       (see removeMean) 
%       minLen   - minimal time series length
%       winSize  - size of the sliding window used to average the local SNR
%
%Output:
%       newTsMap       - processed tsMap (trend, mean, nan removal)
%       Up             - Upper noise limit - cellarray
%       Dw             - Lower noise limit - cellarray
%       snrMap         - snrMap            - cellarray
%       globalSNR      - average SNR for the entire tsMap (scalar)
%       smoothedTsMap  - only used for visualization
%       smoothedSNRMap - only used for visualization
%
% See also: getTimeSeriesNoiseLevel, testImf
%Marco Vilela, 2012

ip=inputParser;
ip.addRequired('tsMap',@(x) isnumeric(x));
ip.addOptional('alpha',   .05,@isscalar);
ip.addOptional('outLevel',7,@isscalar);
ip.addOptional('trend',   1,@isscalar);
ip.addOptional('minLen',  30,@isscalar);
ip.addOptional('winSize', 10,@isscalar);
ip.addOptional('noiseStd', 1,@isscalar);

ip.parse(tsMap,varargin{:});
alpha    = ip.Results.alpha;
outLevel = ip.Results.outLevel;
trend    = ip.Results.trend;
minLen   = ip.Results.minLen;
winSize  = ip.Results.winSize;
noiseStd = ip.Results.noiseStd;

[~,nVar]  = size(tsMap);
%newTsMap  = cell(nVar,1);   
displayTS  = zeros(size(tsMap));
displaySNR = nan(size(tsMap));

for iVar = 1:nVar
    
    newTsMap(iVar) = getTimeSeriesNoiseLevel(tsMap(:,iVar),'alpha',alpha,'winSize',winSize,...
                     'outLevel',outLevel,'trendType',trend,'minLength',minLen,'noiseStd',noiseStd);
    displayTS(newTsMap(iVar).interval,iVar)  = newTsMap(iVar).points;
    displaySNR(newTsMap(iVar).interval,iVar) = newTsMap(iVar).localSNR;

end


%% Smoothing maps
%Parameters for the smoothing operation
smoothParam = .7;
pixSize     = .5;


smoothedTsMap = csaps({1:size(displayTS',1),...
1:size(displayTS',2)},displayTS',smoothParam,{1:pixSize:size(displayTS',1),...
1:pixSize:size(displayTS',2)} );

smoothedSNRMap = csaps({1:size(displaySNR',1),...
1:size(displaySNR',2)},displaySNR',smoothParam,{1:pixSize:size(displaySNR',1),...
1:pixSize:size(displaySNR',2)} );
