function  [procTS,exclude] = timeSeriesPreProcessing(TS,varargin)
%This function performes the following time series operations:
%   - remove outliers
%   - interpolate NaN gaps (be careful with the size of the gap)
%   - remove trend (different types)
%
%Usage:
%       [outTS,exclude] = timeSeriesPreProcessing(TS,varargin)
%
%Input:
%       TS - time series. Format (nVar,nPoints)
%       alpha - alpha level to test signal's imf against white noise imf's. (see getTimeSeriesTrend)
%       nSurr - number of TS surragates created to build the white noise imf distribution
%       minLength - minimal TS length
%       trendType - see getTimeSeriesTrend; -1 for no trend removal
%       gapSize   - length of the NaN gap to be interpolated - usually 1
%       outLevel  - used to detect outliers (see detectOutliers)
%
%Output:
%       procTS   - processed TS
%       exclude - list of variables that didn't pass the minimal length
%                 requirement after processing
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('TS',@(x) ismatrix(x));
[nVar,nObs] = size(TS);

ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('nSurr',100,@isscalar);
ip.addParamValue('minLength',30,@isscalar);
ip.addParamValue('trendType',-1,@isscalar);
ip.addParamValue('gapSize',0,@isscalar);
ip.addParamValue('outLevel',0,@isscalar);
ip.addParamValue('interval',{1:nObs},@iscell);

ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
nSurr    = ip.Results.nSurr;
minLen   = ip.Results.minLength;
trend    = ip.Results.trendType;
gapSize  = ip.Results.gapSize;
interval = ip.Results.interval;
outLevel = ip.Results.outLevel;

procTS   = nan(size(TS));

for iInt = 1:numel(interval)
    
    outTS   = TS(:,interval{iInt});
    exclude = [];
    
    %% Removing trend
    if trend > -1
        
        [auxTS,exclude1] = getTimeSeriesTrend(outTS,'trendType',trend,'nSurr',nSurr,'alpha',alpha,'minLength',minLen);
        outTS            = auxTS.dTS;
        exclude          = union(exclude,exclude1);
    end
    
    %% Removing outliers
    if outLevel > 0
        
        for iVar = 1:nVar
            outTS(iVar,detectOutliers(outTS(iVar,:),outLevel)) = NaN;
        end
        
    end
   
    procTS(:,interval{iInt}) = outTS;
end

%% Interpolating nan Gaps
        %Closing nan gaps <= gapSize .
        %IMPORTANT - Artificial autocorrelation is generated if the gapSize >= 2

if gapSize > 0
    
    for iVar = 1:nVar
        procTS(iVar,:) = gapInterpolation(procTS(iVar,:),gapSize);
    end
    
end

%% Removing by minimum length again in case outliers were detected
exclude1 = find( sum(isfinite(procTS),2) < minLen )';
exclude  = union(exclude1,exclude);