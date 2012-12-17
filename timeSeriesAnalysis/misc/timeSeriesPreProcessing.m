function  [outTS,exclude] = timeSeriesPreProcessing(TS,varargin)
%TS (nVar,nPoints)

ip = inputParser;
ip.addRequired('TS',@(x) isnumeric(x));
ip.addOptional('alpha',.05,@isscalar);
ip.addOptional('nSurr',100,@isscalar);
ip.addOptional('minLength',30,@isscalar);
ip.addOptional('plotYes',0,@isscalar);
ip.addOptional('trendType',0,@isscalar);
ip.addOptional('gapSize',0,@isscalar);
ip.addParamValue('outLevel',7,@isscalar);

ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
nSurr    = ip.Results.nSurr;
minLen   = ip.Results.minLength;
plotYes  = ip.Results.plotYes;
trendT   = ip.Results.trendType;
gapSize  = ip.Results.gapSize;
outLevel = ip.Results.outLevel;


nVar        = size(TS,1);
outTS      = cell(1,nVar);

if outLevel
    TS(detectOutliers(TS,outLevel)) = NaN;
end

for iVar = 1:nVar
    
    %Closing nan gaps <= gapSize . 
    %IMPORTANT - Artificial autocorrelation is generated if the gapSize >= 2
    outTS{iVar} = gapInterpolation(TS(iVar,:),gapSize);
end

%Empty windows 
empIdx  = find( cellfun(@isempty,outTS) );
%Windows with #points < minLength
minIdx  = find( cell2mat( cellfun(@(x) lt(numel(x),minLen),outTS,'UniformOutput',0) ) );
%Windows to be excluded
exclude = unique([empIdx(:);minIdx(:)]);
