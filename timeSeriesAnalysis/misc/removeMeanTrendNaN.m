function [workTS,interval,trend,imf] = removeMeanTrendNaN(TS,varargin)
%Removes mean, trend and NaN from input time series TS
%
%Synopsis:
%         [outTS,interval] = removeMeanTrendNaN(TS)   
%Input:
%       TS        - time series (number of points,number of variables)
%
%       trendType - optional: a scalar giving the type of trend to remove
%                   -1: no trend removal
%                   0 : remove only sample means (see dtrend.m)
%                   1 : default - remove linear trend (see dtrend.m)
%                   2 : remove all deterministic trend
%
%Output:
%       outTS{# of variables}(# of good points)  - cell array with a continuous time series points
%       interval - final interval = initial - (NaN blocks + outliers)  
%       trend    - trend removed from the time series
%       imf      - if trendType is 2, the intrisic mode functions of the corresponding time series
%
%Marco Vilela, 2011

% Input check
ip=inputParser;
ip.addRequired('TS',@isnumeric);
ip.addOptional('trendType',1,@(x)isscalar(x) && ismember(x,-1:2));
ip.parse(TS,varargin{:})
trendType=ip.Results.trendType;

% Initialize output
[nObs,nVar] = size(TS);
workTS      = cell(1,nVar);
interval    = cell(1,nVar);
trend       = cell(1,nVar);
imf         = cell(1,nVar);
count       = 1;
if trendType == 2, imf = cell(1,nVar); end

for i=1:nVar
    
    xi          = find(isnan(TS(:,i)));
    [nanB,nanL] = findBlock(xi,1);
    exclude     = [];
    
%% Excluding gaps larger than 2 points and single NaN extremes
    for j = 1:length(nanB)
        
        if nanL(j) > 2 || ~isempty(intersect(nanB{j},nObs)) || ~isempty(intersect(nanB{j},1))
            
            xi      = setdiff(xi,nanB{j});%Single NaN islands 
            exclude = sort(cat(1,nanB{j},exclude));
            
        end
        
    end

    
    if nObs - length( exclude ) >= 30% forced by the spline used in preWhitening
%% Interpolating single NaN points throughout the time series        
        if ~isempty(xi)%After excluding points, xi is a vector of size 1 NaN 
            
            x         = find(~isnan(TS(:,i)));
            [fB,fL]   = findBlock(union(x,xi));
            [~,idxB]  = max(fL);
            workTS{count} = TS(fB{idxB},i);
            
            workTS{count}(isnan(workTS{count})) = ...
                interp1(intersect(x,fB{idxB}),TS(intersect(x,fB{idxB}),i),intersect(xi,fB{idxB}),'spline');
%% If there are only blocks of NaN(no single NaN), get the largest continuous block of real points            
        else 
            
            [fB,fL]     = findBlock(setdiff(1:nObs,exclude));
            [~,idxB]    = max(fL);
            workTS{count}   = TS(fB{idxB},i);
            
        end
%% Applying the detrend operation after excluding NaN         
        interval{count} = fB{idxB};
        if ismember(trendType,[0 1])
            % Remove sample means or linear trend
            dWorkTS = dtrend(workTS{count},trendType);
            trend{count} = workTS{count} - dWorkTS;
            workTS{count} = dWorkTS;
        elseif trendType == 2
            % Remove deterministic components using preWhitening
            [workTS{count},trend{count},imf{count}]   = preWhitening(workTS{count});
        end
        used(count) = i;
        count = count + 1;
        
    end
    
end

empIdx           = find(cellfun(@isempty,workTS));
workTS(empIdx)   = [];
interval(empIdx) = [];
trend(empIdx)    = [];
imf(empIdx)      = [];