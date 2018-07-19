% [cellMask cellBoundary] = getCellMaskMSS(img, varargin) estimates the cell mask/outline using multi-scale steerable filters
%
% Inputs:
%             img : input image
% 
% Options:
%        'Scales' : vector of scales (sigma) used by the filter. Default: [1 2 4].
%   'FilterOrder' : order of the filters. Default: 3.
%  'RemoveRadius' : radius of the final erosion/refinement step
%
% Outputs:
%        cellMask : binary mask of the cell 
%    cellBoundary : binary mask of the cell outline


% Francois Aguet, September 2011 (last modified: 10/23/2011)


function [cellMask, cellEdge, edgeConfidence] = cellMaskFromMatchedEdges(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Scales', [1 2 4], @isvector);
ip.addParamValue('FilterOrder', 3, @(x) ismember(x, [1 3 5]));
ip.addParamValue('SearchRadius', 6, @isscalar);
ip.addParamValue('NormalizeResponse', false, @islogical);
ip.addParamValue('MaskDilationRadius', 20);
ip.addParamValue('Mask', []);
ip.addParamValue('Display', false, @islogical);
ip.parse(img, varargin{:});
scales = ip.Results.Scales;



cellMask = [];
cellEdge = [];
edgeConfidence = [];

[ny,nx] = size(img);
% image border pixels, in column order CCW
borderIdx = [1:ny 2*ny:ny:(nx-1)*ny nx*ny:-1:(nx-1)*ny+1 (nx-2)*ny+1:-ny:ny+1];
borderMask = zeros(ny,nx);
borderMask(borderIdx) = 1;

fMAD = 1 / norminv(0.75, 0, 1);

%------------------------------------------------------------------------------
% I. Multi-scale steerable filter, edge segments
%------------------------------------------------------------------------------
[res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, ip.Results.FilterOrder, scales);
% figure; imagesc(scaleMap); colormap(gray(256)); axis image; colorbar;
% return
if ip.Results.NormalizeResponse
    res = res ./ filterGauss2D(res, 5);
    % Two alternatives:
    %nms = nonMaximumSuppression(res, theta);
    % -or-
    nms = (nms~=0).*res; % better
end

% Mask of candidate edges
edgeMask = double(bwmorph(nms~=0, 'thin'));
% nms = nms.*edgeMask;

% Break any Y or higher order junctions
% generate list of segments and add associated properties
nn = (imfilter(edgeMask, ones(3), 'same')-1) .* edgeMask;
edgeMask = edgeMask .* (nn<=2);

CC = bwconncomp(edgeMask, 8);

% identify and remove single pixels, update: CC, edgeMask
csize = cellfun(@numel, CC.PixelIdxList);
edgeMask(vertcat(CC.PixelIdxList{csize==1})) = 0;
CC.PixelIdxList(csize==1) = [];
CC.NumObjects = numel(CC.PixelIdxList);
CC.NumPixels = csize(csize>1);

% endpoints of all segments
endpointMask = ((imfilter(edgeMask, ones(3), 'same')-1) .* edgeMask)==1;
pixelIdxList = vertcat(CC.PixelIdxList{:});
endpointIdxList = pixelIdxList(endpointMask(pixelIdxList)==1);
CC.EndpointIdx = mat2cell(endpointIdxList', 1, 2*ones(1,CC.NumObjects));

% computeEdgeSegmentProperties(CC, img, theta);%, varargin)

% sort pixels for each segment from one endpoint to the other
D = bwdistgeodesic(logical(edgeMask), endpointIdxList(1:2:end));
D(isinf(D)) = 0;
pixelOrder = mat2cell(D(pixelIdxList)+1, CC.NumPixels, 1);
[~,pixelOrder] = cellfun(@sort, pixelOrder, 'UniformOutput', false);

for i = 1:CC.NumObjects
    CC.PixelIdxList{i} = CC.PixelIdxList{i}(pixelOrder{i});
end
% CC.avgInt = cellfun(@(px) sum(img(px)), CC.PixelIdxList) ./ CC.NumPixels;
CC.AvgRes = cellfun(@(px) sum(res(px)), CC.PixelIdxList) ./ CC.NumPixels;
% orientation returned by steerable filter is orthogonal to feature
CC.Theta = cellfun(@(px) theta(px)+pi/2, CC.PixelIdxList, 'UniformOutput', false);

% figure; imagesc(edgeMask); colormap(gray(256)); axis image; colorbar;
% hold on;
% plot(
return

% labels of connected components
labels = double(labelmatrix(CC));


%------------------------------------------------------------------------------
% II. Link edge segments based on proximity and angles
%------------------------------------------------------------------------------
%%
[matchedMask] = matchSegmentEndPoints(CC, 'SearchRadius', 4, 'Display', true);
% matchedMask = double(matchedMask);

[matchedMask] = matchSegmentEndPointsOLD(edgeMask, theta, 'SearchRadius', 4, 'Display', true);


%%


% CC = updateEdgeInfo(matchedMask, CC);



return
% compute pixel order and intensity on each side for all detected edge segments
% and get segment endpoint indexes
%%
tic;
CCx = computeEdgeSegmentProperties(CC, img, theta, 'InterpDist1', 1, 'InterpDist2', 10, 'InterpWidth', 2);
toc;
%%
%hval = zeros(1,CC.NumObjects);
%for k = 1:CC.NumObjects
%    hval(k) = kstest2(CC.rval{k}(:), CC.lval{k}(:));
%end

%------------------------------
% Gating of avg int/ avg resp
%------------------------------
tic;
avgRes = cellfun(@(px) sum(res(px)), CC.PixelIdxList) ./ CC.NumPixels;
avgInt = cellfun(@(px) sum(img(px)), CC.PixelIdxList) ./ CC.NumPixels;

[f1,x1] = ksdensity(avgInt, 'npoints', 100);
lmax = locmax1d(f1, 3);
lmin = locmin1d(f1, 3);
idx = find(lmin>lmax(1), 1, 'first');
if ~isempty(idx)
    min0 = lmin(idx);
    T1 = x1(min0);
    T1 = prctile(avgInt(avgInt<T1), 95);
end

[f2,x2] = ksdensity(avgRes(avgInt<T1), 'npoints', 100);
lmax = locmax1d(f2, 3);
lmin = locmin1d(f2, 3);
idx = find(lmin>lmax(1), 1, 'first');
if ~isempty(idx)
    min0 = lmin(idx);
    T2 = x2(min0);
    T2 = prctile(avgRes(avgRes<T2), 95);
else
    T2 = prctile(avgRes, 95);
end
toc;

figure; plot(avgInt, avgRes, 'k.');
hold on;
XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');
plot(T1*[1 1], YLim, 'r');
plot(XLim, T2*[1 1], 'r');

CC2 = updateCC(CC, avgRes>T2 | avgInt>T1);
figure; 
ha(1) = subplot(1,2,1);
%imagesc((labelmatrix(CC2)~=0) .* nms); colormap(gray(256)); axis image; %colorbar;
imagesc((labelmatrix(CC2)~=0)); colormap(gray(256)); axis image; %colorbar
ha(2) = subplot(1,2,2);
imagesc(edgeMask~=0); colormap(gray(256)); axis image; %colorbar;
linkaxes(ha);
return

tmp = zeros(CC.ImageSize); 
for k = 1:CC.NumObjects
    tmp(CC.PixelIdxList{k}) = avgInt(k);
end
figure; imagesc(tmp); colormap(gray(256)); axis image; colorbar;

figure; plot(xi,f);

return
lmax = locmax1d(f, 3);
lmin = locmin1d(f, 3);
% dxi = xi(2)-xi(1);
idx = find(lmin>lmax(1), 1, 'first');
if ~isempty(idx) %&& sum(f(1:lmin(idx(1))))*dxi < modeRatio
    min0 = lmin(idx);
    T = xi(min0);
    CC2 = updateCC(CC, avgIntImg>T);
    figure; imagesc((labelmatrix(CC2)~=0) .* nms); colormap(gray(256)); axis image; colorbar;
    figure; imagesc(nms); colormap(gray(256)); axis image; colorbar;
end



return

%------------------------------------------------------------------------------
% II. Rough estimate of the cell outline based on threshold: coarseMask
%------------------------------------------------------------------------------
% Note: this section is dependent on the robustness of the threshold applied below
coarseMask = ip.Results.Mask;
if isempty(coarseMask)
    % threshold 1st mode (background) of histogram
    img_smooth = filterGauss2D(img, 1);
    T = thresholdFluorescenceImage(img_smooth);
    coarseMask = double(img_smooth>T);
    coarseMask = bwmorph(coarseMask, 'fill'); % clean up isolated negative pixels
end
% get boundary from this mask
bdrVect = bwboundaries(coarseMask);
bdrVect = vertcat(bdrVect{:});
coarseBdr = zeros(ny,nx);
coarseBdr(sub2ind([ny nx], bdrVect(:,1), bdrVect(:,2))) = 1;

% endpoints/intersection of boundary w/ border
borderIS = coarseBdr & borderMask;
borderIS = double(borderIS(borderIdx));
borderIS = borderIdx((conv([borderIS(end) borderIS borderIS(1)], [1 1 1], 'valid')-1)==1);

% clean up, remove image border, add back border intersections
coarseBdr = bwmorph(coarseBdr, 'thin');
coarseBdr(borderIdx) = 0;
coarseBdr(borderIS) = 1;

% generate mask by dilating coarse boundary
coarseBdr = imdilate(coarseBdr, strel('disk', ip.Results.MaskDilationRadius));

% labels within search area
idx = unique(labels.*coarseBdr);
idx(idx==0) = []; % remove background label

% update connected components list
CC.NumObjects = numel(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
CC.NumPixels = CC.NumPixels(idx);
CC.EndpointIdx = CC.EndpointIdx(idx);
CC.pvalProx = CC.pvalProx(idx);
CC.pvalDist = CC.pvalDist(idx);
CC.nvalProx = CC.nvalProx(idx);
CC.nvalDist = CC.nvalDist(idx);

% order interpolated intensities such that the lower intensities are always in 'nval'
for k = 1:CC.NumObjects
    pmean = nanmean(CC.pvalProx{k}(:));
    nmean = nanmean(CC.nvalProx{k}(:));
    if pmean<nmean
        tmp = CC.pvalProx{k};
        CC.pvalProx{k} = CC.nvalProx{k};
        CC.nvalProx{k} = tmp;
        tmp = CC.pvalDist{k};
        CC.pvalDist{k} = CC.nvalDist{k};
        CC.nvalDist{k} = tmp;
    end
end
%------------------------------------------------------------------------------
% Based on coarse boundary, get parameters of background distribution (assumed normal)
M = bwconncomp(~coarseBdr);
mu = cellfun(@(i) mean(img(i)), M.PixelIdxList);
minIdx = find(mu==min(mu), 1, 'first');
mu = mu(minIdx);
sigma = std(img(M.PixelIdxList{minIdx}));

% get robust mean and std of each edge segment
pmean = cellfun(@(i) nanmedian(i(:)), CC.pvalDist);
nmean = cellfun(@(i) nanmedian(i(:)), CC.nvalDist);
pstd = fMAD*cellfun(@(i) mad(i(:),1), CC.pvalDist);
nstd = fMAD*cellfun(@(i) mad(i(:),1), CC.nvalDist);

% average intensity of edge segments
avgIntImg = cellfun(@(px) sum(img(px)), CC.PixelIdxList) ./ CC.NumPixels;

% retain segments with separated intensity distributions on both sides that are not in background
idx = abs(pmean-nmean)>=(pstd+nstd) & (avgIntImg>mu+2*sigma)';
CC = updateCC(CC, idx);


% matrix with average intensity of each segment
avgInt = cellfun(@(px) sum(nms(px)), CC.PixelIdxList) ./ CC.NumPixels;
cellBoundary = zeros(ny,nx);
for k = 1:CC.NumObjects
    cellBoundary(CC.PixelIdxList{k}) = avgInt(k);
end

% Some of the edges are background noise -> bimodal distribution of edge intensities
% % val = nms(edgeMask~=0); % intensities of edges
% % minv = min(val);
% % maxv = max(val);
% % T = graythresh(scaleContrast(val, [], [0 1]));
% % T = T*(maxv-minv)+minv;
% % 

% initial estimate of cell contour
% % cellBoundary = edgeMask > T;
%cellBoundary = edgeMask > min(T, thresholdRosin(val));
%cellBoundary = edgeMask;

% 1st graph matching based on orientation at endpoints, with small search radius
[matchedMask] = matchSegmentEndPoints(cellBoundary, theta,...
    'SearchRadius', 4, 'Display', false);
matchedMask = double(matchedMask);

% merge edge information: interpolated intensity bands on sides, endpoints
CC = updateEdgeInfo(matchedMask, CC);

CC = updateCC(CC, CC.NumPixels>2);
matchedMask = zeros(ny,nx);
matchedMask(vertcat(CC.PixelIdxList{:})) = 1;

% The connected components in this mask are no longer simple segments. For the next matching 
% steps, the two outermost endpoints are needed.
CC.EndpointIdx = getEndpointsGeodesic(matchedMask, horzcat(CC.EndpointIdx{:})');

%-----------------------------
% Identify segments at image border -> candidates for linking to cell edge
%-----------------------------

avgInt = cellfun(@(px) sum(res(px)), CC.PixelIdxList) ./ CC.NumPixels;
CC.AvgInt = avgInt/max(avgInt);
%%
tmp = matchedMask; 
for k = 1:CC.NumObjects
    tmp(CC.PixelIdxList{k}) = CC.AvgInt(k);
end
figure; imagesc(tmp); colormap(gray(256)); axis image; colorbar;


%%
% tmp = matchedMask; tmp([CC.EndpointIdx{:}]) = 2;
% % tmp(vertcat(CC.PixelIdxList{~hval})) = -1;
% figure; imagesc(tmp); colormap(gray(256)); axis image; colorbar;


%----------------------------------------
% Matching based on edge similarity
%----------------------------------------
% TO DO: remove redundant links: only 1 connection per endpoint!
% remove edges that are outside 5th pct of background
matchesFound = true;
borderClosed = false;
iter = 0;
maxiter = 6;
while matchesFound && ~borderClosed && iter<=maxiter;
    [matchLabel, matchIndex] = matchEdgeSegments(CC, 'SearchRadius', 10);
    if isempty(matchIndex)
        matchesFound = false;
    end
    matchedMask = double(linkEdgeSegments(matchedMask, matchIndex));
    CC = updateEdgeInfo(matchedMask, CC);
    avgInt = cellfun(@(px) sum(res(px)), CC.PixelIdxList) ./ CC.NumPixels;
    CC.AvgInt = avgInt/max(avgInt);
    CC.EndpointIdx = getEndpointsGeodesic(matchedMask, horzcat(CC.EndpointIdx{:})');
 
    % check whether endpoints of largest CC are both at image border
    % or whether largest CC is closed contour
    [~,si] = sort(CC.NumPixels, 'descend');
    idx = si(1); % change to more components if several edges in image
    [ye,xe] = ind2sub([ny nx], vertcat(CC.EndpointIdx{idx}));
    if all(xe<=5 | xe>nx-5 | ye<=5 | ye>ny-5)
        borderClosed = true;
    end
    
    iter = iter + 1;
end


% if endpoints close to border, connect
if borderClosed
    [by, bx] = ind2sub([ny,nx], borderIdx);
    [iKD, dist] = KDTreeBallQuery([bx' by'], [xe ye], 5);
    for k = 1:numel(xe)
        if dist{k}(1)~=0 % not self
            [y1, x1] = ind2sub([ny nx], borderIdx(iKD{k}(1)));
            iseg = bresenham([xe(k) ye(k)], [x1 y1]);
            iseg = sub2ind([ny nx], iseg(:,2), iseg(:,1));
            matchedMask(iseg) = 1;
        end
    end
else    
end

% Create new mask with longest edge(s), trim
CC = bwconncomp(matchedMask, 8);
CC = updateCC(CC, idx);
cellEdge = double(labelmatrix(CC));
cellEdge = bwmorph(cellEdge, 'skel');
cellEdge = bwmorph(cellEdge, 'spur', 50);

CC = bwconncomp(~cellEdge, 4);
% sort components by intensity
meanInt = cellfun(@(i) mean(img(i)), CC.PixelIdxList);
[~,si] = sort(meanInt, 'descend');
cellMask = cellEdge;
cellMask(CC.PixelIdxList{si(1)}) = 1;


%----------------------------------------
% Edge confidence
%----------------------------------------
% Boundary coinciding with NMS: 100%
% Closed gaps: response intensity
edgeConfidence = double(cellEdge).*res;
edgeConfidence = edgeConfidence/max(edgeConfidence(:));




if ip.Results.Display
    img0 = scaleContrast(img);
    img1 = img0;
    img0(cellEdge~=0) = 0;
    img1(cellEdge~=0) = 255;
    rgb = uint8(cat(3, img1, img0, img0));
    figure;
    ha(1) = subplot(1,2,1);
    imagesc(img); colormap(gray(256)); axis image;
    ha(2) = subplot(1,2,2);
    imagesc(rgb); colormap(gray(256)); axis image;
    linkaxes(ha);
end






%----------------------------------------
% Matching based on edge strength
%----------------------------------------
% [backMask,backMu,backSig] = estimateBackgroundArea(img);
% bgDist = img(backMask==1);

% check for connected components, whether further matching is required
% tmp = labels==0;
% CCX = bwconncomp(tmp,4)

% % matchesFound = true;
% % iter = 0;
% % %while matchesFound
% %     [matchLabel, matchIndex] = matchEdgeSegments(CC, 'SearchRadius', 10, 'Mode', 'Edge');
% %     if isempty(matchIndex)
% %         matchesFound = false;
% %     end
% %     matchedMask = double(linkEdgeSegments(matchedMask, matchIndex));
% %     CC = updateEdgeInfo(matchedMask, CC);
% %     CC.EndpointIdx = getEndpointsGeodesic(matchedMask, horzcat(CC.EndpointIdx{:})');
% %     %matchLabel
% % 
% %     iter = iter + 1;
% % %end





% update CC based on new mask
function newCC = updateEdgeInfo(matchedMask, CC)

labels = double(labelmatrix(CC));

newCC = bwconncomp(matchedMask, 8);
newCC.NumPixels = cellfun(@numel, newCC.PixelIdxList);
labelMap = mat2cell(labels(vertcat(newCC.PixelIdxList{:})), newCC.NumPixels,1);
labelMap = cellfun(@unique, labelMap, 'UniformOutput', false);

N = newCC.NumObjects;
newCC.EndpointIdx = cell(N,1);
newCC.pvalProx = cell(N,1);
newCC.pvalDist = cell(N,1);
newCC.nvalProx = cell(N,1);
newCC.nvalDist = cell(N,1);
for k = 1:N
    idx = setdiff(labelMap{k}, 0);
    newCC.EndpointIdx{k} = horzcat(CC.EndpointIdx{idx});
    newCC.pvalProx{k} = vertcat(CC.pvalProx{idx});
    newCC.pvalDist{k} = vertcat(CC.pvalDist{idx});
    newCC.nvalProx{k} = vertcat(CC.nvalProx{idx});
    newCC.nvalDist{k} = vertcat(CC.nvalDist{idx});
    %newCC.AvgInt(k) = (CC.AvgInt(idx).*CC.NumPixels(idx))/newCC.NumPixels(k);
end

% update CC based on index
function CC = updateCC(CC, idx)
fnames = fieldnames(CC);
fi = find(structfun(@numel, CC)==CC.NumObjects);
for k = 1:numel(fi)
    CC.(fnames{fi(k)}) = CC.(fnames{fi(k)})(idx);
end
CC.NumObjects = numel(CC.PixelIdxList);

% %------------------------------------------------------------------------------
% % III. Join remaining segments/endpoints using graph matching
% %------------------------------------------------------------------------------
% 
% % Remove long spurs
% cellBoundary = bwmorph(cellBoundary, 'thin');
% cellBoundary = bwmorph(cellBoundary, 'spur', 100);
% cellBoundary = bwmorph(cellBoundary, 'clean'); % spur leaves single pixels -> remove
% 
% % Create mask, use largest connected component within coarse threshold (removes potential loops in boundary)
% maskCC = bwconncomp(~cellBoundary, 4);
% csize = cellfun(@(c) numel(c), maskCC.PixelIdxList);
% [~,idx] = sort(csize, 'descend');
% % two largest components: cell & background
% int1 = mean(img(maskCC.PixelIdxList{idx(1)}));
% int2 = mean(img(maskCC.PixelIdxList{idx(2)}));
% cellMask = zeros(ny,nx);
% if int1 > int2
%     cellMask(maskCC.PixelIdxList{idx(1)}) = 1;
% else
%     cellMask(maskCC.PixelIdxList{idx(2)}) = 1;
% end
% 
% % loop through remaining components, check whether part of foreground or background
% for i = idx(3:end)
%     px = coarseMask(maskCC.PixelIdxList{i});
%     if sum(px) > 0.6*numel(px)
%         cellMask(maskCC.PixelIdxList{i}) = 1;
%     end
% end
% cellMask = imdilate(cellMask, strel('disk',1));
% 
% % Optional: erode filopodia-like structures
% if ~isempty(ip.Results.RemoveRadius)
%     cellMask = imopen(cellMask, strel('disk', ip.Results.RemoveRadius));
% end
%     
% % Final contour: pixels adjacent to mask
% B = bwboundaries(cellMask);
% cellBoundary = zeros(ny,nx);
% cellBoundary(sub2ind([ny nx], B{1}(:,1), B{1}(:,2))) = 1;
% 
% 
% 
% function out = connectEndpoints(inputPoints, queryPoints, radius, labels, cellBoundary, updateBoundary)
% if nargin<6
%     updateBoundary = true;
% end
% 
% dims = size(cellBoundary);
% out = zeros(dims);
% nq = size(queryPoints,1);
% [idx, dist] = KDTreeBallQuery(inputPoints, queryPoints, radius);
% 
% labSelf = labels(sub2ind(dims, queryPoints(:,2), queryPoints(:,1)));
% labAssoc = cellfun(@(i) labels(sub2ind(dims, inputPoints(i,2), inputPoints(i,1))), idx, 'UniformOutput', false);
% 
% % idx of endpoints belonging to other edges
% otherIdx = arrayfun(@(i) labAssoc{i}~=labSelf(i), 1:nq, 'UniformOutput', false);
% 
% % remove segment self-association (and thus query self-association)
% idx = arrayfun(@(i) idx{i}(otherIdx{i}), 1:nq, 'UniformOutput', false);
% dist = arrayfun(@(i) dist{i}(otherIdx{i}), 1:nq, 'UniformOutput', false);
% 
% % generate edge map
% E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:nq, 'UniformOutput', false);
% E = vertcat(E{:});
% 
% if ~isempty(E)
%     idx = E(:,1) < E(:,2);
%     
%     E = E(idx,:); % remove redundancy
%     
%     % generate weights
%     D = vertcat(dist{:});
%     D = D(idx);
%     
%     D = max(D)-D;
%     M = maxWeightedMatching(size(inputPoints,1), E, D);
%     
%     E = E(M,:);
%     
%     % add linear segments corresponding to linked endpoints
%     for i = 1:size(E,1)
%         iseg = bresenham([queryPoints(E(i,1),1) queryPoints(E(i,1),2)],...
%             [inputPoints(E(i,2),1) inputPoints(E(i,2),2)]);
%         out(sub2ind(dims, iseg(:,2), iseg(:,1))) = 1;
%     end
% end
% 
% if updateBoundary
%     out = double(out | cellBoundary);
% end
