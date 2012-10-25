%[mask, matchMask] = matchSegmentEndPoints(mask, theta, varargin) refines a binary edge/ridge map using graph matching to close gaps
%
% Inputs: 
% 
%    mask : binary edge/ridge map (output from steerable filter or similar)
%   theta : orientation map (from steerable filter or similar)
%
% Output:
%
%   matchedMask : refined mask
%  unmatchedIdx : index of unmatched segment endpoints

% Francois Aguet, 01/22/2012 (last modified 10/25/2012)

function [matchedMask, unmatchedIdx] = matchSegmentEndPoints(mask, theta, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('mask');
ip.addRequired('theta');
ip.addParamValue('SearchRadius', 5, @isscalar);
ip.addParamValue('KeepJunctions', true, @islogical);
ip.addParamValue('Display', false, @islogical);
ip.parse(mask, theta, varargin{:});
R = ip.Results.SearchRadius;

mask = double(bwmorph(mask~=0, 'thin'));
[ny, nx] = size(mask);

theta = theta+pi/2;

%=============================================
% I. Generate segments, connected components
%=============================================

% # neighbors
nn = (imfilter(mask, ones(3), 'same')-1) .* mask;

% identify junctions and remove
junctionMatrix = nn>2;
segmentMatrix = mask .* ~junctionMatrix;

% endpoints of all segments
endpointMatrix = (imfilter(segmentMatrix, ones(3), 'same')-1) .* segmentMatrix;
endpointMatrix = endpointMatrix==1;

% generate list of segments and add associated properties
CC = bwconncomp(segmentMatrix, 8);

% identify and remove single pixels, update segment matrix
issingleton = cellfun(@(i) numel(i)==1, CC.PixelIdxList);
CC.NumObjects = CC.NumObjects - sum(issingleton);
singletonIdx = CC.PixelIdxList(issingleton);
segmentMatrix([singletonIdx{:}]) = 0;
CC.PixelIdxList(issingleton) = [];

% labels of connected components
labels = double(labelmatrix(CC));

ns = CC.NumObjects; % # segments

% add endpoint indexes
endpointIdx = find(endpointMatrix);
for i = 1:ns
    CC.endpointIdx{i} = intersect(CC.PixelIdxList{i}, endpointIdx);
    %CC.isSegment(i) = max(nn(CC.PixelIdxList{i}))<3;
end
% sorted index
endpointIdx = [CC.endpointIdx{:}];
endpointIdx = endpointIdx(:);

% order pixels from one endpoint to the other
% CC = computeSegmentProperties(CC, theta);


%=============================================
% II. Match segments
%=============================================
[ye, xe] = ind2sub([ny nx], endpointIdx);

matchedMask = zeros(ny,nx);
matchesFound = true;

% iterative matching
unmatchedIdx = 1:numel(endpointIdx);
matchedIdx = [];

iter = 1;
while matchesFound
    
    % segment endpoints admissible for matching
    X = [xe(unmatchedIdx) ye(unmatchedIdx)];
    [idx, dist] = KDTreeBallQuery(X, X, R);
    
    % Generate all possible pairs from the points that resulted from the query.
    % These pairs are the edges in the graph used for matching
    E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:numel(unmatchedIdx), 'UniformOutput', false);
    E = vertcat(E{:});
    
    if ~isempty(idx)
        
        % remove redundant pairs
        E = E(E(:,1) < E(:,2),:);
        
        % Omitting this in the first pass adds stability by essentially ignoring
        % (self matching) small segments adjacent to (and offset from) gaps
        if iter>1
            % remove queries on same segment
            Elabel = ceil(unmatchedIdx(E)/2);
            E(Elabel(:,1)==Elabel(:,2),:) = [];
        end
        
        % perform matching
        % weights based on angle
        t1 = theta(endpointIdx(unmatchedIdx(E(:,1))));
        t2 = theta(endpointIdx(unmatchedIdx(E(:,2))));
        a1 = abs(t1 - t2);
        a2 = abs(a1-pi);
        minAngle = min(a1,a2);
        cost = cos(minAngle);
        
        % endpoint vectors
        v1 = [cos(t1) sin(t1)]';
        v2 = [cos(t2) sin(t2)]';
        % angle btw. vectors
        dt = acos(sum(v1.*v2,1));
        v2(:,dt>pi/2) = -v2(:,dt>pi/2); % flip if wrong direction
        % vectors between endpoint pairs
        vL = [X(E(:,2),1)-X(E(:,1),1) X(E(:,2),2)-X(E(:,1),2)]';
        vL = vL./repmat(sqrt(sum(vL.^2,1)), [2 1]);
        % mean btw. endpoint vectors
        vMean = (v1+v2) ./ repmat(sqrt(sum((v1+v2).^2,1)), [2 1]);
        diffT = acos(sum(vMean.*vL,1));
        diffT(diffT>pi/2) = pi-diffT(diffT>pi/2);
        
        % The angle difference could be included in the matching cost.
        % Currently, any pair with diffT>pi/4 is discarded
        rmIdx = diffT>pi/4;
        E(rmIdx,:) = [];
        cost(rmIdx) = [];
        
        M = maxWeightedMatching(numel(unmatchedIdx), E, cost); % returns index (M==true) of matches
        % retain only pairs that are matches
        E = E(M,:);
        
        % remove 'self' matches (only relevant for 1st iteration)
        E(labels(endpointIdx(unmatchedIdx(E(:,1))))==labels(endpointIdx(unmatchedIdx(E(:,2)))), :) = [];
        
        if ip.Results.Display
            figure; imagesc(segmentMatrix- 4*matchedMask); colormap(gray(256)); axis image; colorbar;
            hold on;
            plot(xe(unmatchedIdx), ye(unmatchedIdx), 'rx');
            T = theta(endpointIdx(unmatchedIdx));%+pi/2; % +pi/2 only for ridges -> debug
            quiver(X(:,1), X(:,2), cos(T), sin(T),0);
            %quiver(X(E(:,1),1), X(E(:,1),2), vL(1,:)', vL(2,:)',0, 'c');
            
            plot(X(unique(E(:)),1), X(unique(E(:)),2), 'go');
            plot([X(E(:,1),1) X(E(:,2),1)]', [X(E(:,1),2) X(E(:,2),2)]', 'y-')
        end
        
        % fill mask
        for i = 1:size(E,1)
            iseg = bresenham([X(E(i,1),1) X(E(i,1),2)], [X(E(i,2),1) X(E(i,2),2)]);
            matchedMask(sub2ind([ny nx], iseg(:,2), iseg(:,1))) = 1;
        end
        
        matchedIdx = unmatchedIdx(unique(E(:))); % current iter only
        unmatchedIdx = setdiff(unmatchedIdx, matchedIdx);
    end
    
    iter = iter + 1;
    if size(E,1)==0
        matchesFound = false;
    end
end

matchedMask = segmentMatrix | matchedMask;
if ip.Results.KeepJunctions
    matchedMask = matchedMask | junctionMatrix;
end

unmatchedIdx = endpointIdx(unmatchedIdx);

if ip.Results.Display
    figure; imagesc(2*matchedMask-segmentMatrix); colormap(gray(256)); axis image; colorbar;
end
