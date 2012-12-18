function  [outTS,exclude] = timeSeriesPreProcessing(TS,varargin)
%TS (nVar,nPoints)

ip = inputParser;
ip.addRequired('TS',@(x) ismatrix(x));
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('nSurr',100,@isscalar);
ip.addParamValue('minLength',30,@isscalar);
ip.addParamValue('plotYes',0,@isscalar);
ip.addParamValue('trendType',0,@isscalar);
ip.addParamValue('gapSize',0,@isscalar);
ip.addParamValue('outLevel',7,@isscalar);
ip.addParamValue('missingObs',0,@(x) ge(x,0) & le(x,100));

ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
nSurr    = ip.Results.nSurr;
minLen   = ip.Results.minLength;
plotYes  = ip.Results.plotYes;
trendT   = ip.Results.trendType;
gapSize  = ip.Results.gapSize;
outLevel = ip.Results.outLevel;
missObs  = ip.Results.missingObs;

[nVar,nObs] = size(TS);
outTS       = cell(1,nVar);

if outLevel > 0
    TS(detectOutliers(TS,outLevel)) = NaN;
end

for iVar = 1:nVar
    %Closing nan gaps <= gapSize . 
    %IMPORTANT - Artificial autocorrelation is generated if the gapSize >= 2
    outTS{iVar} = gapInterpolation(TS(iVar,:),gapSize);
end

excludeFun = @(x) le(numel(isfinite(x)),minLen*(100 - missObs)/100);
%Windows with #points < minLength
exclude = find( cell2mat( cellfun(@(x) excludeFun(x),outTS,'UniformOutput',0) ) );
