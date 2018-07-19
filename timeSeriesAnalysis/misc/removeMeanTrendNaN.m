function [workTS,interval,trend,imf,excludeVar] = removeMeanTrendNaN(TS,varargin)
%Removes mean, trend and NaN from input time series TS
%
%Synopsis:
%         [outTS,interval,trend,imf,excludeVar] = removeMeanTrendNaN(TS)   
%Input:
%       TS        - time series (number of points,number of variables)
%
%       trendType - optional: a scalar giving the type of trend to remove
%                   -1: no trend removal
%                   0 : remove only sample means (see dtrend.m)
%                   1 : default - remove linear trend (see dtrend.m)
%                   2 : remove all deterministic trend
%
%       minLength  - minimal length accepted. Any window that has less than
%                    minLength will be discarded.
%
%Output:
%       outTS{# of variables}(# of good points)  - cell array with a continuous time series points
%       interval   - final interval = initial - (NaN blocks + outliers)  
%       trend      - trend removed from the time series
%       imf        - if trendType is 2, the intrisic mode functions of the corresponding time series
%       excludeVar - variables that did not pass the minimal length test
%
%See also: preWhitening
%Marco Vilela, 2011

 
% Input check
ip=inputParser;
ip.addRequired('TS',@isnumeric);
ip.addParamValue('trendType',1,@(x)isscalar(x) && ismember(x,-1:2));
ip.addParamValue('minLength',30,@(x)isscalar(x));

ip.parse(TS,varargin{:})
trendType = ip.Results.trendType;
minLength = ip.Results.minLength;

% Initialize output
nVar        = size(TS,2);
workTS      = cell(1,nVar);
interval    = cell(1,nVar);
trend       = cell(1,nVar);
imf         = cell(1,nVar);
exC         = 1;
excludeVar  = [];

 
if trendType == 2, imf = cell(1,nVar); end

 
for iVar = 1:nVar
    
%% Interpolating single NaN points throughout the time series            
     workTS{iVar}  = gapInterpolation(TS(:,iVar),1);
    [xAxisB,bLeng] = findBlock(find(isfinite(workTS{iVar})),1);
    [~,idxB]       = max(bLeng);    
    
    if bLeng(idxB) >= minLength 
        
        % Applying the detrend operation after excluding NaN         
        interval{iVar} = xAxisB{idxB};
        if ismember(trendType,[0 1])
            
            % Remove sample means or linear trend
            dWorkTS      = detrend(workTS{iVar}(interval{iVar}),trendType);
            trend{iVar}  = workTS{iVar}(interval{iVar}) - dWorkTS;
            workTS{iVar} = dWorkTS;
            
        elseif trendType == 2
            
            % Remove deterministic components using preWhitening
            [workTS{iVar},trend{iVar},imf{iVar}] = preWhitening(workTS{iVar}(interval{iVar}));
        else
            workTS{iVar} = workTS{iVar}(interval{iVar});
        end
        
    else
        
        excludeVar(exC) = iVar;
        exC = exC + 1;
        
    end

    
end

%Empty windows 
empIdx           = find( cellfun(@isempty,workTS) );
%Windows with #points < minLength
minIdx           = find( cell2mat( cellfun(@(x) lt(numel(x),minLength),workTS,'UniformOutput',0) ) );
%All windows to be excluded
finIdx           = union(empIdx,minIdx);

excludeVar       = union(excludeVar,finIdx);

workTS(excludeVar) = [];
interval(excludeVar) = [];
trend(excludeVar) = [];
imf(excludeVar) = [];