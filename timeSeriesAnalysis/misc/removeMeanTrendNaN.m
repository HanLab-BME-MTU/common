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
ip.addOptional('trendType',1,@(x)isscalar(x) && ismember(x,0:2));
ip.parse(TS,varargin{:})
trendType=ip.Results.trendType;

% Initialize output
[nobs,nvar] = size(TS);
workTS      = cell(1,nvar);
interval    = cell(1,nvar);
trend       = cell(1,nvar);
if trendType == 2, imf = cell(1,nvar); end

for i=1:nvar
    
    xi          = find(isnan(TS(:,i)));
    [nanB,nanL] = findBlock(xi,1);
    exclude     = [];
    
    for j = 1:length(nanB)%excluding gaps larger than 2 points and extremes (1 and N points)
        
        if nanL(j) > 2 || ~isempty(intersect(nanB{j},nobs)) || ~isempty(intersect(nanB{j},1))
            
            xi      = setdiff(xi,nanB{j});
            exclude = sort(cat(1,nanB{j},exclude));
            
        end
        
    end
    
    if ~isempty(xi)%After excluding points, xi is a vector of 1 NaN block
        
        x         = find(~isnan(TS(:,i)));
        [fB,fL]   = findBlock(union(x,xi));
        [~,idxB]  = max(fL);
        workTS{i} = TS(fB{idxB},i);
        
        workTS{i}(isnan(workTS{i})) = ...
            interp1(intersect(x,fB{idxB}),TS(intersect(x,fB{idxB}),i),intersect(xi,fB{idxB}),'spline');
        
    elseif (nobs - length(exclude) ) >= 4 % forced by the spline used in preWhitening
        
        [fB,fL]     = findBlock(setdiff(1:nobs,exclude));
        [~,idxB]    = max(fL);
        workTS{i}   = TS(fB{idxB},i);
        
    end
    
    interval{i} = fB{idxB};
    if ismember(trendType,[0 1])
        % Remove sample means or linear trend
        dWorkTS = dtrend(workTS{i},trendType);
        trend{i} = workTS{i} - dWorkTS;
        workTS{i} = dWorkTS;
    elseif trendType==2
        % Remove deterministic components using preWhitening
        [workTS{i},trend{i},imf{i}]   = preWhitening(workTS{i});
    end
end