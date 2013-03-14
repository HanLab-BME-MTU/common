function  [outTS,exclude] = timeSeriesPreProcessing(TS,varargin)
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
%       trendType - see getTimeSeriesTrend
%       gapSize   - length of the NaN gap to be interpolated - usually 1
%       outLevel  - used to detect outliers (see detectOutliers)
%
%Output:
%       outTS   - resulting TS after operations
%       exclude - list of variables that were excluded because some of the minimal length requirement
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('TS',@(x) ismatrix(x));
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('nSurr',100,@isscalar);
ip.addParamValue('minLength',30,@isscalar);
ip.addParamValue('trendType',-1,@isscalar);
ip.addParamValue('gapSize',0,@isscalar);
ip.addParamValue('outLevel',0,@isscalar);

ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
nSurr    = ip.Results.nSurr;
minLen   = ip.Results.minLength;
trendT   = ip.Results.trendType;
gapSize  = ip.Results.gapSize;
outLevel = ip.Results.outLevel;


nVar  = size(TS,1);
outTS = TS;
%% Removing trend
if trendT > -1
    
    auxTS = getTimeSeriesTrend(outTS,'trendType',trendT,'nSurr',nSurr,'alpha',alpha);
    outTS = auxTS.dTS;   
    
end


%% Removing outliers
if outLevel > 0
    
    for iVar = 1:nVar
        outTS(iVar,detectOutliers(TS(iVar,:),outLevel)) = NaN;
    end
    
end

%% Interpolating nan Gaps
        %Closing nan gaps <= gapSize .
        %IMPORTANT - Artificial autocorrelation is generated if the gapSize >= 2

if gapSize > 0
    
    for iVar = 1:nVar
        outTS(iVar,:) = gapInterpolation(TS(iVar,:),gapSize);
    end
    
end

%% Removing by minimum length
exclude = find( sum(isfinite(outTS),2) < minLen );