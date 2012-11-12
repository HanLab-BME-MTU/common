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
ip.addOptional('trendType',1,@(x)isscalar(x) && ismember(x,-1:2));
ip.addOptional('minLength',30,@(x)isscalar(x));

ip.parse(TS,varargin{:})
trendType = ip.Results.trendType;
minLength = ip.Results.minLength;

% Initialize output
[nObs,nVar] = size(TS);
workTS      = cell(1,nVar);
interval    = cell(1,nVar);
trend       = cell(1,nVar);
imf         = cell(1,nVar);
count       = 1;
exC         = 1;
excludeVar  = [];

 
if trendType == 2, imf = cell(1,nVar); end

 
for i = 1:nVar

    
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
    
%% Interpolating single NaN points throughout the time series            
    if nObs - length( exclude ) >= minLength% forced by the spline used in preWhitening

        if ~isempty(xi)%After excluding points, xi is a vector of size 1 NaN 
          
            x         = find(~isnan(TS(:,i)));
            [fB,fL]   = findBlock(union(x,xi));
            [~,idxB]  = max(fL);
            workTS{count} = TS(fB{idxB},i);

            
            workTS{count}(isnan(workTS{count})) = ...
                interp1(intersect(x,fB{idxB}),TS(intersect(x,fB{idxB}),i),intersect(xi,fB{idxB}),'spline');
            
        else % If there are only blocks of NaN(no single NaN), get the largest continuous block of real points            

            [fB,fL]     = findBlock(setdiff(1:nObs,exclude));
            [~,idxB]    = max(fL);
            workTS{count}   = TS(fB{idxB},i);

        end
        
        % Applying the detrend operation after excluding NaN         
        interval{count} = fB{idxB};
        if ismember(trendType,[0 1])
            
            % Remove sample means or linear trend
            dWorkTS       = dtrend(workTS{count},trendType);
            trend{count}  = workTS{count} - dWorkTS;
            workTS{count} = dWorkTS;
            
        elseif trendType == 2
            
            % Remove deterministic components using preWhitening
            [workTS{count},trend{count},imf{count}]   = preWhitening(workTS{count});
            
        end
        
        count = count + 1;
        
    else
        
        excludeVar(exC) = i;
        exC = exC + 1;
        
    end

    
end

 
empIdx           = find( cellfun(@isempty,workTS) );
minIdx           = find( cell2mat( cellfun(@(x) lt(numel(x),minLength),workTS,'UniformOutput',0) ) );
finIdx           = union(empIdx,minIdx);
workTS(finIdx)   = [];
interval(finIdx) = [];
trend(finIdx)    = [];
imf(finIdx)      = [];
excludeVar       = union(excludeVar,finIdx);