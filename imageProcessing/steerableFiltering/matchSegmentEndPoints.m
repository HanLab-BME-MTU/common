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

% Francois Aguet, 01/22/2012 (last modified 10/24/2012)

function [matchedMask] = matchSegmentEndPoints(mask, theta, varargin)

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
end
% sorted index
endpointIdx = [CC.endpointIdx{:}];
endpointIdx = endpointIdx(:);


% order pixels from one endpoint to the other
[ny, nx] = size(mask);
getNH = @(i) [i-ny-1 i-ny i-ny+1 i-1 i+1 i+ny-1 i+ny i+ny+1]; 

for i = 1:ns
    np = numel(CC.PixelIdxList{i});
    unorderedIdx = CC.PixelIdxList{i};
    orderedIdx = zeros(1,np);
    
    % assign 1st endpoint (choice is arbitrary)
    orderedIdx(1) = CC.endpointIdx{i}(1);
    unorderedIdx(unorderedIdx==orderedIdx(1)) = [];
    
    for k = 2:np
        % local neighborhood indexes
        hoodIdx = getNH(orderedIdx(k-1));
        nextIdx = intersect(unorderedIdx, hoodIdx);
        unorderedIdx(unorderedIdx==nextIdx) = [];
        orderedIdx(k) = nextIdx;
    end
    CC.PixelIdxList{i} = orderedIdx;
    CC.rawAngle{i} = theta(orderedIdx);
    % running average - POTENTIAL BUG: angle wrapping in (-pi/2..pi/2]
%     if np>w
%         CC.smoothAngle{i} = conv([CC.rawAngle{i}(b:-1:2) CC.rawAngle{i} CC.rawAngle{i}(end-1:-1:end-b+1)], ones(1,w)/w, 'valid');
%     else
%         CC.smoothAngle{i} = CC.rawAngle{i};
%     end
end

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
    
    iter = iter + 1;
    if size(E,1)==0
        matchesFound = false;
    end
end

matchedMask = segmentMatrix | matchedMask;
if ip.Results.KeepJunctions
    matchedMask = matchedMask | junctionMatrix;
end

if ip.Results.Display
    figure; imagesc(2*matchedMask-segmentMatrix); colormap(gray(256)); axis image; colorbar;
end


% junk, erase in future version

% [ny nx] = size(mask);
% dims = [ny nx];
% borderIdx = [1:ny (nx-1)*ny+(1:ny) ny+1:ny:(nx-2)*ny+1 2*ny:ny:(nx-1)*ny];
% 
% 
% % minimal connectivity
% mask = double(bwmorph(mask~=0, 'thin'));
% 
% % endpoints
% nn = padarrayXT(mask, [1 1]);
% nn = conv2(sumKernel, sumKernel', nn, 'valid');
% nn = (nn-1).*mask;
% [ye,xe] = ind2sub([ny nx], find(nn==1));
% 
% % connected components
% CC = bwconncomp(mask, 8);
% labels = double(labelmatrix(CC));
% 
% % admissible endpoints
% [idx, dist] = KDTreeBallQuery([xe ye], [xe ye], ip.Results.MaxRadius);
% 
% labSelf = labels(sub2ind([ny nx], ye, xe));
% labAssoc = cellfun(@(i) labels(sub2ind([ny nx], ye(i), xe(i))), idx, 'UniformOutput', false);
% 
% % idx of endpoints belonging to other edges
% % otherIdx = arrayfun(@(i) labAssoc{i}~=labSelf(i), 1:numel(xe), 'UniformOutput', false);
% % 
% % % remove segment self-association (and thus query self-association)
% % idx = arrayfun(@(i) idx{i}(otherIdx{i}), 1:numel(xe), 'UniformOutput', false);
% % dist = arrayfun(@(i) dist{i}(otherIdx{i}), 1:numel(xe), 'UniformOutput', false);
% 
% 
% 
% 
% % generate edge map
% E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:numel(xe), 'UniformOutput', false);
% E = vertcat(E{:});
% 
% idx = E(:,1) < E(:,2);
% 
% E = E(idx,:); % remove redundancy
% 
% % weights based on angle
% T = theta(sub2ind([ny nx], ye(E(:,1)), xe(E(:,1)))) - theta(sub2ind([ny nx], ye(E(:,2)), xe(E(:,2))));
% T = abs(cos(T));
% 
% 
% % weights based on distance
% % D = vertcat(dist{:});
% % D = D(idx);
% % D = max(D)-D;
% M = maxWeightedMatching(numel(xe), E, T);
% E = E(M,:);
% 
% % remove matches from same segment
% l1 = labels(sub2ind([ny nx], ye(E(:,1)), xe(E(:,1))));
% l2 = labels(sub2ind([ny nx], ye(E(:,2)), xe(E(:,2))));
% E(l1==l2,:) = [];
% 
% matchMask = zeros(ny,nx);
% 
% for i = 1:size(E,1)
%     iseg = bresenham([xe(E(i,1)) ye(E(i,1))], [xe(E(i,2)) ye(E(i,2))]);
%     matchMask(sub2ind(dims, iseg(:,2), iseg(:,1))) = 1;
% end
% 
% mask = double(mask | matchMask);
% 
% mask = bwmorph(mask, 'close');
% mask = double(bwmorph(mask, 'thin'));
% mask(borderIdx) = 1;
% mask = double(bwmorph(mask, 'spur', 100));
% mask = bwmorph(mask, 'clean'); % spur leaves single pixels -> remove
% mask(borderIdx) = 0;





% update labels
% CC = bwconncomp(mask, 8);
% labels = double(labelmatrix(CC));

% dilate each endpoint with 3x3 square SE
% epMask = zeros(dims);
% epMask(sub2ind(dims, ye, xe)) = 1;
% % dilMask = imdilate(epMask, strel('square', 3));
% dilMask = bwmorph(epMask, 'dilate');
