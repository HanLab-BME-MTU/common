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

function [cellMask cellBoundary] = getCellMaskMSS(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Scales', [1 2 4], @isvector);
ip.addParamValue('FilterOrder', 3, @(x) ismember(x, [1 3 5]));
ip.addParamValue('RemoveRadius', [], @isscalar);
ip.addParamValue('Mask', []);
ip.parse(img, varargin{:});
scales = ip.Results.Scales;

[ny,nx] = size(img);
borderIdx = [1:ny (nx-1)*ny+(1:ny) ny+1:ny:(nx-2)*ny+1 2*ny:ny:(nx-1)*ny];

%------------------------------------------------------------------------------
% I. Multi-scale steerable filter
%------------------------------------------------------------------------------
ns = length(scales);
res = cell(1,ns);
th = cell(1,ns);
nms = cell(1,ns);
for si = 1:ns
    [res{si}, th{si}, nms{si}] = steerableDetector(img, ip.Results.FilterOrder, scales(si));
end
%maxIdx = ones(ny,nx);
maxRes = res{1};
maxNMS = nms{1};
%maxTh = th{1};
for si = 2:ns
    idx = res{si}>maxRes;
    maxRes(idx) = res{si}(idx);
    %maxIdx(idx) = si;
    maxNMS(idx) = nms{si}(idx);
    %maxTh(idx) = th{si}(idx);
end

% Mask of candidate edges
maxNMS = maxNMS.*bwmorph(maxNMS~=0, 'thin'); % or skel

% Break any Y or higher order junctions
nn = padarrayXT(double(maxNMS~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
nn = (nn-1) .* (maxNMS~=0);
junctionMask = nn>2;
maxNMS(junctionMask) = 0;

% Individual edges: connected components
CC = bwconncomp(maxNMS,8);
csize = cellfun(@(c) numel(c), CC.PixelIdxList);

% Remove singletons
nsingle = sum(csize==1);
CC.NumObjects = CC.NumObjects-nsingle;
CC.PixelIdxList(csize==1) = [];
csize(csize==1) = [];


%------------------------------------------------------------------------------
% II. Rough estimate of the cell outline based on threshold: bdrMask
%------------------------------------------------------------------------------
mask = ip.Results.Mask;
if isempty(mask)
    % threshold 1st mode (background) of histogram
    img_smooth = filterGauss2D(img, 1);
    T = thresholdFluorescenceImage(img_smooth);
    mask = double(img_smooth>T);
    mask = bwmorph(mask, 'fill');
end
bdrVect = bwboundaries(mask);
bdrMask = zeros(ny,nx);
% largest perimeter/boundary
plength = cellfun(@(i) size(i,1), bdrVect);
idx = plength==max(plength);
bdrMask(sub2ind([ny nx], bdrVect{idx}(:,1), bdrVect{idx}(:,2))) = 1;
bdrMask = bwmorph(bdrMask, 'thin');
bdrMask = bwmorph(bdrMask, 'spur', 1);

% Remove pixels at image border
bdrMask(borderIdx) = 0;

% Area of admissible edges based on bdrMask; labels of these edges
bdrMask = imdilate(bdrMask, strel('disk', 20)); % dilation is arbitrary...
idx = unique(double(labelmatrix(CC)).*bdrMask);
idx(idx==0) = [];
CC.NumObjects = numel(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);
csize = csize(idx);

% Mask of retained edges
edgeMask = labelmatrix(CC) ~= 0;

% average intensity of each connected component
avgInt = cellfun(@(px) sum(maxNMS(px)), CC.PixelIdxList) ./ csize;
avgMap = zeros(ny,nx);
for k = 1:CC.NumObjects
    avgMap(CC.PixelIdxList{k}) = avgInt(k);
end

% Some of the edges are background noise -> bimodal distribution of edge intensities
val = maxNMS(edgeMask); % intensities of edges
minv = min(val);
maxv = max(val);
T = graythresh(scaleContrast(val, [], [0 1]));
T = T*(maxv-minv)+minv;
% initial estimate of cell contour
cellBoundary = avgMap > T;

edgeLabels = double(labelmatrix(CC));

for i = 1:2
    % junctions connected to edges above threshold
    cjunctions = bwmorph(cellBoundary, 'dilate') .* junctionMask;
    
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
    px = mask(maskCC.PixelIdxList{i});
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
        %iseg = bresenham([xe(E(i,1)) ye(E(i,1))], [xe(E(i,2)) ye(E(i,2))]);
        iseg = bresenham([queryPoints(E(i,1),1) queryPoints(E(i,1),2)],...
            [inputPoints(E(i,2),1) inputPoints(E(i,2),2)]);
        out(sub2ind(dims, iseg(:,2), iseg(:,1))) = 1;
    end
end

if updateBoundary
    out = double(out | cellBoundary);
end
