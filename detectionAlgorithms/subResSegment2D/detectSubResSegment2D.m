function [params, ImBG] = detectSubResSegment2D(ima,mask,sigmaPSF,minSize,bitDepth)

%% --- Initialization --- %%

% Make sure image's class is double
if ~isa(ima,'double')
    ima = double(ima);
end

[nrows,ncols] = size(ima);

%% --- Step 1 ---- %%

% Filter image with 2nd order steerable filter
% TODO: Use M=4 instead
filteredIma1 = steerableFiltering(ima,2,sigmaPSF);
filteredIma1(filteredIma1 < 0) = 0;
filteredIma1(mask == false) = 0;

% Filter image with 2-D laplacian filter
filteredIma2 = filterLoG(ima,sigmaPSF);
filteredIma2(mask == false) = 0;

%% --- Step 2 ---- %%

% Segment filteredIma1
% TODO: User LSM
segmentMask = logical(blobSegmentThreshold(filteredIma1,minSize,0,mask));

segmentLabels = bwlabel(segmentMask);
CCstats = regionprops(segmentMask,'PixelIdxList','BoundingBox');
nCC = numel(CCstats);

%% --- Step 3 ---- %%

% Get local maxima of filteredIma2 that belong to BW
% Note: we could have set hside to the full size of an individual speckle
% (2 * ceil(3 * sigma) + 1) but it is too restrictive since some speckles
% can be closer and somehow overlap with each other.
hside = 2 * ceil(sigmaPSF) + 1;
locMax = locmax2d(filteredIma2,[hside hside]);
locMax(segmentMask == false) = 0;
indPSF = find(locMax ~= 0);
nPSF = numel(indPSF);

%% --- Step 4 ---- %%

cands(1:numel(indPSF)) = struct('Lmax',[],'IBkg',[],'status',[]);
for iPSF = 1:nPSF
    [y x] = ind2sub([nrows ncols], indPSF(iPSF));
    cands(iPSF).Lmax = [y x];
    cands(iPSF).IBkg = mean(ima(CCstats(segmentLabels(indPSF(iPSF))).PixelIdxList)) / (2^bitDepth-1);
    cands(iPSF).status = true;
end

stdNoise = std(ima(segmentMask == false & mask == true) / (2^bitDepth-1));

%alpha.alphaR = .05;
alpha.alphaA = .01;
% alpha.alphaD = .1;
alpha.alphaF = 0;

% uncomment
%finalPoints = detectSubResFeatures2D(ima,cands,sigmaPSF,alpha,0,1,bitDepth,0,stdNoise);

%% --- Step 5 ---- %%
tic;
% Possible optimization: instead of computed Radon transform of
% filteredIma1 for every possible angle, use the Feature-adapted Radon
% transform.

winSize = 2*ceil(2*sigmaPSF)+1;

thetaRadon = 0:179;
theta = -thetaRadon*pi/180;
ct = cos(theta);
st = sin(theta);

for iCC = 1:nCC
    % Floor bounding box
    bb = ceil(CCstats(iCC).BoundingBox);

    % Get the footprint of the CC
    CC = zeros(bb(4),bb(3));
    [y x] = ind2sub([nrows ncols],CCstats(iCC).PixelIdxList);
    x = x - bb(1) + 1;
    y = y - bb(2) + 1;
    indLocal = sub2ind(size(CC), y, x);
    CC(indLocal) = 1;
    
    % Compute the Radon transform on the CC's footprint. It gives the
    % length of each Radon lines.
    [L xp] = radon(CC, thetaRadon);
    
    % Crop R using bounding box and restrinct to CC's footprint
    filteredCC1 = zeros(bb(4),bb(3));
    filteredCC1(indLocal) = filteredIma1(CCstats(iCC).PixelIdxList);
    
    % Compute the Radon transform on the crop filtered image
    RT = radon(filteredCC1, thetaRadon);
    
    % Compute the mean integral of R along the radon lines.
    averageRT = RT ./ (L + 1e-10);
    
    % The Radon origin is floor((size(BW)+1)/2).
    cRadon = floor((bb(3:4)+1)/2);
    
    % Get the local maxima per column (for each orientation) of averageRT
    locMax = locmax2d(averageRT,[winSize 1]);
    locMax(averageRT <= 0) = 0;
    locMax = locMax ~= 0;
    
    %% Compute initial value of the center's coordinates of each segment candidates
    
    % Compute the distance of each local maxima away from the main axis,
    % for each main axis's orientation
    distAside = bsxfun(@times,locMax,xp);
    
    % Compute the distance of each local maxima away form the secondary
    % axis, for each secondary axis's orientation. For that, compute the
    % Radon transform of a function defines as the signed distance from the
    % perpendical axis passing through cRadon. The orientation of that
    % perpendical follows each radon orientation.
    D1 = zeros(size(CC));
    D1(indLocal) = cRadon(1) - x;
    D2 = zeros(size(CC));
    D2(indLocal) = y - cRadon(2);
    RD = bsxfun(@times, radon(D1, thetaRadon), st) + ...
        bsxfun(@times, radon(D2, thetaRadon), ct);
    
    distAlong = zeros(size(locMax));
    distAlong(locMax) = RD(locMax) ./ L(locMax);
    
    segmentX = bsxfun(@times,distAside,ct) + bsxfun(@times,distAlong,-st);
    segmentX(locMax) = bb(1) + cRadon(1) - 1 + segmentX(locMax);
    segmentY = bsxfun(@times,distAside,st) + bsxfun(@times,distAlong, ct);
    segmentY(locMax) = bb(2) + cRadon(2) - 1 + segmentY(locMax);
    
    %% Compute initial value of the amplitude of each segment %%
    
    % Crop ima using bounding box and restrinct to CC's footprint
    imaCC = zeros(bb(4),bb(3));
    imaCC(indLocal) = ima(CCstats(iCC).PixelIdxList);
    
    % Compute the Radon transform on the crop image. averageRI correspond
    % to the mean intensity of the image along the radon lines.
    RI = radon(imaCC, thetaRadon);
    
    segmentAmp = zeros(size(locMax));
    segmentAmp(locMax) = RI(locMax) ./ L(locMax);
    
    %% Compute initial value of the length of each segment candidates
    
    segmentLength = zeros(size(locMax));
    segmentLength(locMax) = L(locMax);
    
    %% Sequential Parameter Optimization
    
    for iTheta = 1:numel(theta)
        % For each orientation, optimize the segment parameter one after
        % the other
    
        % Keep orientation that minimize the norm of residual
    end
    
    %% Optimization of all parameters
    
    % After the optimization, test the residual using alphaR quantile and
    % plot segment in 2 colors whether they pass the test or not.
end
toc
%% --- Step 6 --- %%

%% --- Step 7 --- %%

%% --- Step 8 --- %%

params = [];
