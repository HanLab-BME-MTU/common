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
filteredIma2 = filterLaplacian2D(ima,sigmaPSF);
filteredIma2(mask == false) = 0;

%% --- Step 2 ---- %%

% Segment filteredIma1
% TODO: User LSM
BW = logical(blobSegmentThreshold(filteredIma1,minSize,0,mask));

%% --- Step 3 ---- %%

% Get local maxima of filteredIma2 that belong to BW
% Note: we could have set hside to the full size of an individual speckle
% (2 * ceil(3 * sigma) + 1) but it is too restrictive since some speckles
% can be closer and somehow overlap with each other.
hside = 2 * ceil(sigmaPSF) + 1;
locMax = locmax2d(filteredIma2,[hside hside]);
locMax(BW == false) = 0;
indPSF = find(locMax ~= 0);
nPSF = numel(indPSF);

%% --- Step 4 ---- %%

L = bwlabel(BW);
CCstats = regionprops(BW,'PixelIdxList');

cands(1:numel(indPSF)) = struct('Lmax',[],'IBkg',[],'status',[]);
for iPSF = 1:nPSF
    [y x] = ind2sub([nrows ncols], indPSF(iPSF));
    cands(iPSF).Lmax = [y x];
    cands(iPSF).IBkg = min(ima(CCstats(L(indPSF(iPSF))).PixelIdxList)) / (2^bitDepth-1);
    cands(iPSF).status = true;
end

stdNoise = std(ima(BW == false & mask == true) / (2^bitDepth-1));

finalPoints = detectSubResFeatures2D(ima,cands,sigmaPSF,[],0,0,bitDepth,0,stdNoise);

%% --- Step 5 ---- %%


%% --- Step 6 ---- %%


%% --- Step 7 ---- %%


%% --- Step 8 ---- %%
