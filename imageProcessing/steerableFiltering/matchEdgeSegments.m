% Francois Aguet, 10/2012

function [matchLabel, matchIndex] = matchEdgeSegments(CC, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('SearchRadius', 10, @isscalar);
ip.addParamValue('Mode', 'KSDistance', @(x) any(strcmpi(x, {'KSDistance', 'Edge'})));
ip.parse(varargin{:});

dims = CC.ImageSize;
labels = double(labelmatrix(CC));

pidx = vertcat(CC.PixelIdxList{:});

endpointIdx = vertcat(CC.EndpointIdx{:});
endpointLabels = labels(endpointIdx);

[yi, xi] = ind2sub(dims, pidx);
X = [xi yi];
[yi, xi] = ind2sub(dims, endpointIdx);
[idx, dist] = KDTreeBallQuery(X, [xi yi], ip.Results.SearchRadius);

np = numel(endpointIdx);
matchLabel = NaN(np,2);
matchIndex = NaN(np,2);
matchDist = NaN(np,1);
% loop through endpoints, determine matches
for k = 1:np
    % remove self-queries
    rmIdx = labels(pidx(idx{k}))==endpointLabels(k);
    idx{k}(rmIdx) = [];
    dist{k}(rmIdx) = [];
    % label of closest point
    if ~isempty(idx{k})
        iMatchLabel = [endpointLabels(k) labels(pidx(idx{k}(1)))];
        iMatchIndex = [endpointIdx(k) pidx(idx{k}(1))];
        [~,sidx] = sort(iMatchLabel);
        matchLabel(k,:) = iMatchLabel(sidx);
        matchIndex(k,:) = iMatchIndex(sidx);
        matchDist(k) = dist{k}(1);
    end
end
% determine unique pairs (shortest distance)
[~,sidx] = sort(matchDist);
matchLabel = matchLabel(sidx,:);
matchIndex = matchIndex(sidx,:);
rmIdx = isnan(matchLabel(:,1));
matchLabel(rmIdx,:) = [];
matchIndex(rmIdx,:) = [];
[~,sidx] = sort(matchLabel(:,1));
matchLabel = matchLabel(sidx,:);
matchIndex = matchIndex(sidx,:);

N = size(matchLabel,1);
valid = false(N,1);
valid(1) = true;
for k = 2:N
    valid(k) = ~any(matchLabel(1:k-1,1)==matchLabel(k,1) & matchLabel(1:k-1,2)==matchLabel(k,2));
end
matchLabel = matchLabel(valid,:);
matchIndex = matchIndex(valid,:);
N = size(matchLabel,1);

switch ip.Results.Mode
    case 'KSDistance'
        % cost based on KS distance
        cost = zeros(N,1);
        for k = 1:N
            [~,~,ksLL] = kstest2(CC.lval{matchLabel(k,1)}(:), CC.lval{matchLabel(k,2)}(:));
            [~,~,ksHH] = kstest2(CC.rval{matchLabel(k,1)}(:), CC.rval{matchLabel(k,2)}(:));
            
            %[~,~,ks1] = kstest2(CC.lval{matchList(k,1)}(:), CC.rval{matchList(k,1)}(:));
            %[~,~,ks2] = kstest2(CC.lval{matchList(k,2)}(:), CC.rval{matchList(k,2)}(:));
            
            cost(k) = 1-max([ksLL ksHH]);
            %cost(k) = 1-mean([ksLL ksHH]);
            
            % penalize cost when left/right distributions of a segment are close
            %cost(k) = cost(k) * min(ks1,ks2);
            %cost(k) = cost(k) * (1-(1-ks1)*(1-ks2));
            
        end
        rmIdx = cost<0.2;
        matchLabel(rmIdx,:) = [];
        matchIndex(rmIdx,:) = [];
        cost(rmIdx) = [];
    case 'Edge'
        % score edges: background vs. foreground
        medDiff = cellfun(@(i) nanmedian(i(:)), CC.rval) - cellfun(@(i) nanmedian(i(:)), CC.lval);
        medDiff = (medDiff-min(medDiff)) / (max(medDiff)-min(medDiff));
        nsize = CC.NumPixels/sum(CC.NumPixels); % ! geodesic distance should be used here !
        % cost = mean(medDiff(matchList),2)
        % reward strong and long edges
        cost = sum(nsize(matchLabel),2).*mean(medDiff(matchLabel),2);
        %rmIdx = cost<0.1;
        % matchList(rmIdx,:) = [];
        % cost(rmIdx) = [];
end

M = maxWeightedMatching(CC.NumObjects, matchLabel, cost); % returns index (M==true) of matches
matchLabel = matchLabel(M,:);
matchIndex = matchIndex(M,:);
