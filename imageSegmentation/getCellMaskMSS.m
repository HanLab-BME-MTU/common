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

function [cellMask, cellBoundary] = getCellMaskMSS(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Scales', [1 2 4], @isvector);
ip.addParamValue('FilterOrder', 3, @(x) ismember(x, [1 3 5]));
ip.addParamValue('SearchRadius', 6, @isscalar);
ip.addParamValue('NormalizeResponse', false, @islogical);
ip.addParamValue('Mask', []);
ip.parse(img, varargin{:});
scales = ip.Results.Scales;

[ny,nx] = size(img);
% ordered index, column-order CCW
borderIdx = [1:ny 2*ny:ny:(nx-1)*ny nx*ny:-1:(nx-1)*ny+1 (nx-2)*ny+1:-ny:ny+1];
borderMask = zeros(ny,nx);
borderMask(borderIdx) = 1;

%------------------------------------------------------------------------------
% I. Multi-scale steerable filter
%------------------------------------------------------------------------------
[res, theta, nms] = multiscaleSteerableDetector(img, ip.Results.FilterOrder, scales);

if ip.Results.NormalizeResponse
    res = res ./ filterGauss2D(res, 5);
    %nms = nonMaximumSuppression(res, theta);
    % -or-
    nms = (nms~=0).*res;
end

% Mask of candidate edges
edgeMask = double(bwmorph(nms~=0, 'thin'));

% Break any Y or higher order junctions
nn = (imfilter(edgeMask, ones(3), 'same')-1) .* edgeMask;

junctionMatrix = nn>2;
segmentMatrix = edgeMask .* ~junctionMatrix;

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
csize = cellfun(@numel, CC.PixelIdxList);

% labels of connected components
labels = double(labelmatrix(CC));

%------------------------------------------------------------------------------
% II. Rough estimate of the cell outline based on threshold: coarseMask
%------------------------------------------------------------------------------
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
borderIS = borderIdx((conv([borderIS(end) borderIS borderIS(1)], [1 1 1], 'same')-1)==1);

% clean up, remove image border, add intersects
coarseBdr = bwmorph(coarseBdr, 'thin');
coarseBdr(borderIdx) = 0;
coarseBdr(borderIS) = 1;

edgeSearchMask = imdilate(coarseBdr, strel('disk', 20)); % dilation is arbitrary...

% labels within search area
idx = unique(labels.*edgeSearchMask);
idx(idx==0) = []; % remove background label

% update connected components list
CC.NumObjects = numel(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
csize = csize(idx);

% mask with average intensity of each segment
avgInt = cellfun(@(px) sum(nms(px)), CC.PixelIdxList) ./ csize;
edgeMask = zeros(ny,nx);
for k = 1:CC.NumObjects
    edgeMask(CC.PixelIdxList{k}) = avgInt(k);
end

% Some of the edges are background noise -> bimodal distribution of edge intensities
val = nms(edgeMask~=0); % intensities of edges
minv = min(val);
maxv = max(val);
T = graythresh(scaleContrast(val, [], [0 1]));
T = T*(maxv-minv)+minv;

% initial estimate of cell contour
cellBoundary = edgeMask > T;
% cellBoundary = edgeMask;

% 1st graph matching based on orientation at endpoints, with small search radius
[matchedMask] = matchSegmentEndPoints(cellBoundary, theta, 'SearchRadius', ip.Results.SearchRadius, 'Display', false);

% find endpoint candidates on skeleton
nn = double(bwmorph(matchedMask, 'thin'));
nn = (imfilter(nn, ones(3), 'same')-1) .* nn;
endpointMatrix = nn==1;
endpointIdx = find(endpointMatrix);

CC = bwconncomp(matchedMask, 8);
for k = 1:CC.NumObjects
    CC.endpointIdx{k} = intersect(CC.PixelIdxList{k}, endpointIdx);
end

% retain the two endpoints that are furthest apart for later matching
nEndpoint = cellfun(@numel, CC.endpointIdx);
for k = 1:max(nEndpoint)
    idx = find(nEndpoint>=max(3,k));
    % seed point indexes
    sp = cellfun(@(x) x(k), CC.endpointIdx(idx));
    seedMatrix = false(ny,nx);
    seedMatrix(sp) = true;
    D = bwdistgeodesic(matchedMask, seedMatrix);
    for i = 1:numel(idx)
        CC.endpointDist{idx(i)}(k) = max(D(CC.PixelIdxList{idx(i)}));
    end    
end
idx = find(nEndpoint>2);
for k = 1:numel(idx)
    [~,si] = sort(CC.endpointDist{idx(k)}, 'descend');
    [~,si] = sort(si, 'descend');
    CC.endpointIdx{idx(k)} = CC.endpointIdx{idx(k)}(si<=2);
end
CC = rmfield(CC, 'endpointDist');

matchedMask = double(matchedMask);

csize = cellfun(@numel, CC.PixelIdxList);
avgInt = cellfun(@(px) sum(res(px)), CC.PixelIdxList) ./ csize;
for k = 1:CC.NumObjects
    matchedMask(CC.PixelIdxList{k}) = avgInt(k);
    %CC.isSegment(k) = max(nn(CC.PixelIdxList{k}))<3;
    %CC.endpointIdx{k} = intersect(CC.PixelIdxList{k}, endpointIdx);
    if isempty(CC.endpointIdx{k})
        CC.endpointIdx{k} = CC.PixelIdxList{k}(1);
    end
end

CC = computeSegmentProperties(CC, img, theta);

figure; imagesc(rgbOverlay(img, matchedMask, [1 0 0])); colormap(gray(256)); axis image; colorbar;

cellMask = [];
return






edgeLabels = double(labelmatrix(CC));
for i = 1:2
    % junctions connected to edges above threshold
    cjunctions = bwmorph(cellBoundary, 'dilate') .* junctionMatrix;
    
    % add segments connected to these junctions
    clabels = unique(double(bwmorph(cjunctions, 'dilate')) .* edgeLabels);
    cellBoundary = cellBoundary | cjunctions | ismember(edgeLabels, clabels(2:end));
end
cellBoundary = double(cellBoundary);

% Endpoints of edges
endpoints = double((cellBoundary .* (conv2(sumKernel, sumKernel', padarrayXT(cellBoundary, [1 1]), 'valid')-1))==1);

% last step of hysteresis: recover edges close to these endpoints
labels = unique(bwmorph(endpoints, 'dilate', 2) .* edgeLabels);
cellBoundary = cellBoundary | ismember(edgeLabels, labels(2:end));

% extend endpoints that are within 1 pixel of image border
border = zeros(ny,nx);
border(borderIdx) = 1;
endpoints(borderIdx) = 0;
cellBoundary = cellBoundary | (bwmorph(endpoints, 'dilate') & border);
cellBoundary = double(bwmorph(cellBoundary~=0, 'thin'));

% update endpoints
endpoints = double((cellBoundary .* (conv2(sumKernel, sumKernel', padarrayXT(cellBoundary, [1 1]), 'valid')-1))==1);

%------------------------------------------------------------------------------
% III. Join remaining segments/endpoints using graph matching
%------------------------------------------------------------------------------
CC = bwconncomp(cellBoundary, 8);
labels = double(labelmatrix(CC));

% A) Link endpoints that are in close proximity
[ye,xe] = find(endpoints~=0);
cellBoundary = connectEndpoints([xe ye], [xe ye], 1.42*2, labels, cellBoundary);

% B) Some junctions are not detected by the filters -> join endpoints with nearby edges
endpoints = double((cellBoundary .* (conv2(sumKernel, sumKernel', padarrayXT(cellBoundary, [1 1]), 'valid')-1))==1);
CC = bwconncomp(cellBoundary, 8);
labels = double(labelmatrix(CC));
[ye,xe] = find(endpoints~=0);
[yb,xb] = find(cellBoundary~=0);

cellBoundary = connectEndpoints([xb yb], [xe ye], 1.42*2, labels, cellBoundary);

% C) For each remaining connected component, identify the most distant (geodesic distance) endpoints
CC = bwconncomp(cellBoundary, 8);
labels = double(labelmatrix(CC));

epLabels = zeros(ny,nx);
for c = 1:CC.NumObjects
    % endpoints for this segment
    endpoints = double((cellBoundary .* (conv2(sumKernel, sumKernel', padarrayXT(double(labels==c), [1 1]), 'valid')-1))==1);
    endpoints(borderIdx) = 0;
    [yi,xi] = find(endpoints~=0);
    nep = numel(xi);
    if nep > 2
        % use each endpoint as a seed point
        xm = zeros(1,nep);
        ym = zeros(1,nep);
        dm = zeros(1,nep);
        for i = 1:nep
            D = bwdistgeodesic(labels==c, xi(i), yi(i));
            dm(i) = max(D(:));
            [ym(i) xm(i)] = ind2sub([ny nx], find(D==dm(i), 1', 'first'));
        end
        idx = find(dm==max(dm), 1, 'first');
        epLabels(yi(idx), xi(idx)) = c;
        epLabels(ym(idx), xm(idx)) = c;
    else
        epLabels(sub2ind([ny nx], yi, xi)) = c;
    end
end

% D) Bridge longer gaps linearly (real gap only if no intersections)
[ye,xe] = find(epLabels~=0);
segments = connectEndpoints([xe ye], [xe ye], 20, epLabels, cellBoundary, false);

segCC = bwconncomp(segments);
% read pixel positions in boundary, count. If count >2, intersection
validSeg = cellfun(@(i) sum(cellBoundary(i)), segCC.PixelIdxList) == 2;

% update boundary
cellBoundary(vertcat(segCC.PixelIdxList{validSeg})) = 1;


% if a single segment is left, but the endpoints are not at the image boundary, join them
epLabels(borderIdx) = 0;
nep = sum(epLabels(:)~=0);
if CC.NumObjects==1 && nep == 2
    [yi,xi] = find(epLabels~=0);
    seg = bresenham([xi(1) yi(1)], [xi(2) yi(2)]);
    cellBoundary(sub2ind([ny nx], seg(:,2), seg(:,1))) = 1;
end


% Remove long spurs
cellBoundary = bwmorph(cellBoundary, 'thin');
cellBoundary = bwmorph(cellBoundary, 'spur', 100);
cellBoundary = bwmorph(cellBoundary, 'clean'); % spur leaves single pixels -> remove

% Create mask, use largest connected component within coarse threshold (removes potential loops in boundary)
maskCC = bwconncomp(~cellBoundary, 4);
csize = cellfun(@(c) numel(c), maskCC.PixelIdxList);
[~,idx] = sort(csize, 'descend');
% two largest components: cell & background
int1 = mean(img(maskCC.PixelIdxList{idx(1)}));
int2 = mean(img(maskCC.PixelIdxList{idx(2)}));
cellMask = zeros(ny,nx);
if int1 > int2
    cellMask(maskCC.PixelIdxList{idx(1)}) = 1;
else
    cellMask(maskCC.PixelIdxList{idx(2)}) = 1;
end

% loop through remaining components, check whether part of foreground or background
for i = idx(3:end)
    px = coarseMask(maskCC.PixelIdxList{i});
    if sum(px) > 0.6*numel(px)
        cellMask(maskCC.PixelIdxList{i}) = 1;
    end
end
cellMask = imdilate(cellMask, strel('disk',1));

% Optional: erode filopodia-like structures
if ~isempty(ip.Results.RemoveRadius)
    cellMask = imopen(cellMask, strel('disk', ip.Results.RemoveRadius));
end
    
% Final contour: pixels adjacent to mask
B = bwboundaries(cellMask);
cellBoundary = zeros(ny,nx);
cellBoundary(sub2ind([ny nx], B{1}(:,1), B{1}(:,2))) = 1;



function out = connectEndpoints(inputPoints, queryPoints, radius, labels, cellBoundary, updateBoundary)
if nargin<6
    updateBoundary = true;
end

dims = size(cellBoundary);
out = zeros(dims);
nq = size(queryPoints,1);
[idx, dist] = KDTreeBallQuery(inputPoints, queryPoints, radius);

labSelf = labels(sub2ind(dims, queryPoints(:,2), queryPoints(:,1)));
labAssoc = cellfun(@(i) labels(sub2ind(dims, inputPoints(i,2), inputPoints(i,1))), idx, 'UniformOutput', false);

% idx of endpoints belonging to other edges
otherIdx = arrayfun(@(i) labAssoc{i}~=labSelf(i), 1:nq, 'UniformOutput', false);

% remove segment self-association (and thus query self-association)
idx = arrayfun(@(i) idx{i}(otherIdx{i}), 1:nq, 'UniformOutput', false);
dist = arrayfun(@(i) dist{i}(otherIdx{i}), 1:nq, 'UniformOutput', false);

% generate edge map
E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:nq, 'UniformOutput', false);
E = vertcat(E{:});

if ~isempty(E)
    idx = E(:,1) < E(:,2);
    
    E = E(idx,:); % remove redundancy
    
    % generate weights
    D = vertcat(dist{:});
    D = D(idx);
    
    D = max(D)-D;
    M = maxWeightedMatching(size(inputPoints,1), E, D);
    
    E = E(M,:);
    
    % add linear segments corresponding to linked endpoints
    for i = 1:size(E,1)
        iseg = bresenham([queryPoints(E(i,1),1) queryPoints(E(i,1),2)],...
            [inputPoints(E(i,2),1) inputPoints(E(i,2),2)]);
        out(sub2ind(dims, iseg(:,2), iseg(:,1))) = 1;
    end
end

if updateBoundary
    out = double(out | cellBoundary);
end
