function  [outTS,exclude] = timeSeriesPreProcessing(TS,varargin)
%This function performe the following time series operations:
%   - remove outliers
%   - interpolate NaN gaps (be careful with the size of the gap)
%   - remove trend (different types)
%
%TS (nVar,nPoints)

ip = inputParser;
ip.addRequired('TS',@(x) ismatrix(x));
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('nSurr',100,@isscalar);
ip.addParamValue('minLength',30,@isscalar);
ip.addParamValue('trendType',-1,@isscalar);
ip.addParamValue('gapSize',0,@isscalar);
ip.addParamValue('outLevel',7,@isscalar);

ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
nSurr    = ip.Results.nSurr;
minLen   = ip.Results.minLength;
trendT   = ip.Results.trendType;
gapSize  = ip.Results.gapSize;
outLevel = ip.Results.outLevel;


nVar  = size(TS,1);
outTS = TS;

%% Removing outliers
if outLevel > 0
    
    outTS(detectOutliers(TS,outLevel)) = NaN;
    
end

%% Interpolating nan Gaps
if gapSize > 0
    
    for iVar = 1:nVar
        %Closing nan gaps <= gapSize .
        %IMPORTANT - Artificial autocorrelation is generated if the gapSize >= 2
        outTS(iVar,:) = gapInterpolation(TS(iVar,:),gapSize);
    end
    
end

%% Removing trend
if trendType > -1
    
    outTS = getTimeSeriesTrend(outTS,'trendType',trendT,'nSurr',nSurr,'alpha',alpha);
        
end

%% Removing by minimum length
exclude = find( sum(isfinite(outTS),2) < minLen );