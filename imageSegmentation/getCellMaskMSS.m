% [cellMask cellBoundary] = getCellMaskMSS(img, varargin) estimates the cell mask/outline using multi-scale steerable filters
%
% Inputs:
%             img : input image
% 
% Options:
%        'Scales' : vector of scales (sigma) used by the filter. Default: [1 2 4].
%   'FilterOrder' : order of the filters. Default: 3.
%
% Outputs:
%        cellMask : binary mask of the cell 
%    cellBoundary : binary mask of the cell outline

% Francois Aguet, September 2011 (last modified: 09/28/2011)

function [cellMask cellBoundary] = getCellMaskMSS(img, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addParamValue('Scales', [1 2 4], @isvector);
ip.addParamValue('FilterOrder', 3, @(x) ismember(x, [1 3 5]));
ip.addParamValue('Mask', []);
ip.parse(img, varargin{:});
scales = ip.Results.Scales;

[ny,nx] = size(img);

% Multi-scale steerable filter
ns = length(scales);
res = cell(1,ns);
th = cell(1,ns);
nms = cell(1,ns);
for si = 1:ns
    [res{si}, th{si}, nms{si}] = steerableDetector(img, ip.Results.FilterOrder, scales(si));
end
maxIdx = ones(ny,nx);
maxRes = res{1};
maxNMS = nms{1};
maxTh = th{1};
for si = 2:ns
    idx = res{si}>maxRes;
    maxRes(idx) = res{si}(idx);
    maxIdx(idx) = si;
    maxNMS(idx) = nms{si}(idx);
    maxTh(idx) = th{si}(idx);
end

% Mask of candidate edges
nmsMask = bwmorph(maxNMS~=0, 'thin'); % skel
nmsTest = maxNMS.*nmsMask;

% Break any Y or higher order junctions
nn = padarrayXT(double(nmsMask), [1 1]);
k = [1 1 1];
nn = conv2(k, k', nn, 'valid');
nn = (nn-1).*nmsMask;
nmsMask(nn>2) = 0;
nmsTest = nmsTest.*nmsMask;


% Individual edges: connected components
CC = bwconncomp(nmsMask,8);
csize = cellfun(@(c) numel(c), CC.PixelIdxList);

% Remove singletons
nsingle = sum(csize==1);
CC.NumObjects = CC.NumObjects-nsingle;
CC.PixelIdxList(csize==1) = [];
csize(csize==1) = [];

avgInt = cellfun(@(px) sum(nmsTest(px)), CC.PixelIdxList) ./ csize;

% generate corresponding images
avgMap = zeros(ny,nx);
for k = 1:CC.NumObjects
    avgMap(CC.PixelIdxList{k}) = avgInt(k);
end



% Initial estimate of the cell outline
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
bdrMask(sub2ind([ny nx], bdrVect{1}(:,1), bdrVect{1}(:,2))) = 1;
bdrMask = bwmorph(bdrMask, 'thin');
bdrMask = bwmorph(bdrMask, 'spur', 1);

% Remove pixels at image border
bdrMask(:,[1 end]) = 0;
bdrMask([1 end],:) = 0;

% Dilate to get mask of suitable edges; get labels of these edges
bdrMask = imdilate(bdrMask, strel('disk', 10)); % dilation is arbitrary...
idx = unique(double(labelmatrix(CC)).*bdrMask);
idx(idx==0) = [];
CC.NumObjects = numel(idx);
CC.PixelIdxList = CC.PixelIdxList(idx);

% Mask of retained edges
fmask = labelmatrix(CC) ~= 0;

% Some of the edges are background noise -> bimodal distribution of edge intensities
val = maxRes(fmask); % intensities of edges
minv = min(val);
maxv = max(val);
T = graythresh(scaleContrast(val, [], [0 1]));
T = T*(maxv-minv)+minv;
finalCand = fmask.*avgMap;
finalCand(finalCand<T) = 0; 

% High-confidence edges
labelsFinal = double(labelmatrix(CC));
labelsFinal = unique(labelsFinal(finalCand~=0));
CC.NumObjects = numel(labelsFinal);
CC.PixelIdxList = CC.PixelIdxList(labelsFinal);

% The previous threshold excludes some weaker edge segments. These are rescued successively
% using morphological operations, thus avoiding a hysteresis threshold-based approach.

% Two levels of 'rescue': junctions, then connected objects
junctionMask = nn>2;
junctionCC = bwconncomp(junctionMask,8); 
validLabels = unique(double(bwmorph(finalCand, 'dilate')) .* double(labelmatrix(junctionCC)));
validLabels(validLabels==0) = [];
junctionCC.NumObjects = numel(validLabels);
junctionCC.PixelIdxList = junctionCC.PixelIdxList(validLabels);
cellBoundary = finalCand>0 | double(labelmatrix(junctionCC));

% Retrieve connected objects from NMS before threshold T
CC = bwconncomp(fmask.*avgMap | labelmatrix(junctionCC), 8);
labelsFinal = unique(double(labelmatrix(CC)) .* double(bwmorph(cellBoundary, 'dilate', 2)));
labelsFinal(labelsFinal==0) = [];
CC.NumObjects = numel(labelsFinal);
CC.PixelIdxList = CC.PixelIdxList(labelsFinal);
cellBoundary = double(labelmatrix(CC)>0);

% Endpoints of edges
k = [1 1 1];
endpoints = double((cellBoundary .* (conv2(k, k', padarrayXT(cellBoundary, [1 1]), 'valid')-1))==1);
endpoints(:,[1 end]) = 0;
endpoints([1 end],:) = 0;

% Dilate endpoints; keep points with at least two neighbors with different labels
tmp = double(bwmorph(endpoints, 'dilate')) - endpoints;
tmp(:,[1 end]) = 0;
tmp([1 end],:) = 0;

[yi,xi] = find(tmp~=0);
ne = numel(xi);
labels = double(labelmatrix(CC));
nlabels = zeros(1,ne);
for k = 1:ne
    nlabels(k) = numel(unique(labels(yi(k)-1:yi(k)+1,xi(k)-1:xi(k)+1)));% includes 0!
end
idx = nlabels > 2;
links = zeros(ny,nx);
links(sub2ind([ny nx], yi(idx),xi(idx))) = 1;

% Reduce each component to 1 pixel
linkCC = bwconncomp(links,8);
k = [1 1 1];
nn = (conv2(k, k', padarrayXT(links, [1 1]), 'valid')-1).*links;
for k = 1:linkCC.NumObjects
    nni = nn(linkCC.PixelIdxList{k});
    linkCC.PixelIdxList{k} = linkCC.PixelIdxList{k}(nni==max(nni));
end
cellBoundary = cellBoundary | labelmatrix(linkCC)>0;

% Remove long spurs
cellBoundary = bwmorph(cellBoundary, 'spur', 100);
cellBoundary = bwmorph(cellBoundary, 'clean');

% Create mask, use largest connected component within coarse threshold
maskCC = bwconncomp(~cellBoundary, 4);
csize = cellfun(@(c) numel(c), maskCC.PixelIdxList);
[~,idx] = sort(csize, 'descend');
idx = idx(1:2); % indexes of two largest components
int1 = mean(img(maskCC.PixelIdxList{idx(1)}));
int2 = mean(img(maskCC.PixelIdxList{idx(2)}));
cellMask = zeros(ny,nx);
if int1 > int2
    cellMask(maskCC.PixelIdxList{idx(1)}) = 1;
else
    cellMask(maskCC.PixelIdxList{idx(2)}) = 1;
end
cellMask = imdilate(cellMask, strel('disk',1));

% Final contour: pixels adjacent to mask
cellBoundary = (cellMask+cellBoundary)==2;
