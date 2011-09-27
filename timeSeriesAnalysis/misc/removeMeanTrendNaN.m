function [outTS,interval]=removeMeanTrendNaN(TS,trendType)
%Removes mean, trend and NaN from time series
%Input:
%       TS        - time series (number of points,number of variables)
%       trendType - 'c' to remove mean
%                   'l' to remove linear trend;Default = 'c'
%Output:
%       outTS{# of variables}(# of good points)  - cell array with a continuous time series points
%       interval - good points interval    
%Marco Vilela, 2011

if nargin < 2
    trendType ='c';
end

[nobs,nvar] = size(TS);
workTS      = cell(1,nvar);
interval    = cell(1,nvar);
%***********************
for i=1:nvar
    xi          = find(isnan(TS(:,i)));
    [nanB,nanL] = findBlock(xi,1);
    exclude     = [];
    for j=1:length(nanB)
        if nanL(j) > 2 || ~isempty(intersect(nanB{j},nobs)) || ~isempty(intersect(nanB{j},1))
            xi      = setdiff(xi,nanB{j});
            exclude = sort(cat(1,nanB{j},exclude));
        end
        
    end
    if ~isempty(xi)
        x         = find(~isnan(TS(:,i)));
        [fB,fL]   = findBlock(union(x,xi));
        [~,idxB]  = max(fL);
        workTS{i} = TS(fB{idxB},i);
        workTS{i}(isnan(workTS{i})) = ...
            interp1(intersect(x,fB{idxB}),TS(intersect(x,fB{idxB}),i),intersect(xi,fB{idxB}),'spline');
    elseif length(exclude) < nobs
        [fB,fL]   = findBlock(setdiff(1:nobs,exclude));
        [~,idxB]  = max(fL);
        workTS{i} = TS(fB{idxB},i);
    end
    if length(exclude) < nobs
        interval{i} = fB{idxB};
        workTS{i}   = workTS{i} - repmat(mean(workTS{i}),sum(~isnan(workTS{i})),1);
        workTS{i}   = detrend(workTS{i},trendType);
    end
end
outTS = workTS(cellfun(@(x) ~isempty(x),workTS))';
