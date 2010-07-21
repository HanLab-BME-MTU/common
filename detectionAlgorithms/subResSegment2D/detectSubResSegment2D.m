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

    % Expand bounding box so that every optimized segment model can fit
    % into the expanded bounding box.
    offset = 3 * sigmaPSF;
    for i = 1:2
        bb(i) = min(bb(i) - offset, 1);
        if bb(i) + bb(i+2) + 2 * offset < ncols
            bb(i+2) = bb(i+2) + 2 * offset;
        else
            bb(i+2) = ncols - bb(i) + 1;
        end
    end
    
    % Get the footprint of the CC
    CC = zeros(bb(4),bb(3));
    [y x] = ind2sub([nrows ncols],CCstats(iCC).PixelIdxList);
    x = x - bb(1) + 1 + offset;
    y = y - bb(2) + 1 + offset;
    indLocal = sub2ind(size(CC), y, x);
    CC(indLocal) = 1;
    
    % Compute the Radon transform on the CC's footprint. It gives the
    % length of each Radon lines.
    [L xp] = radon(CC, thetaRadon);
    
    % Crop ima using bounding box and restrinct to CC's footprint
    imaCC = zeros(bb(4),bb(3));
    imaCC(indLocal) = ima(CCstats(iCC).PixelIdxList);
    
    % Compute the Radon transform on the crop image. averageRI correspond
    % to the mean intensity of the image along the radon lines.
    RI = radon(imaCC, thetaRadon);
    averageRI = RI ./ (L + 1e-10);
    
    % Crop R using bounding box and restrinct to CC's footprint
    filteredCC1 = zeros(bb(4),bb(3));
    filteredCC1(indLocal) = filteredIma1(CCstats(iCC).PixelIdxList);
    
    % Compute the Radon transform on the crop filtered image
    RT = radon(filteredCC1, thetaRadon);
    
    % Compute the mean integral of R along the radon lines.
    averageRT = RT ./ (L + 1e-10);
    
    % The Radon origin is floor((size(BW)+1)/2).
    cRadon = floor((bb(3:4)+1)/2);
        
    % Compute the Radon transform of a function defines as the signed
    % distance from the perpendical axis passing through cRadon. The
    % orientation of that perpendical follows each radon orientation.
    D1 = zeros(size(CC));
    D1(indLocal) = cRadon(1) - x;
    D2 = zeros(size(CC));
    D2(indLocal) = y - cRadon(2);
    RD = bsxfun(@times, radon(D1, thetaRadon), st) + bsxfun(@times, radon(D2, thetaRadon), ct);
    
    locMax = locmax2d(averageRT,[winSize 1]);
    locMax(averageRT <= 0) = 0;
    
    % TODO: try to get rid of that loop (is it faster?)
    for iTheta = 1:numel(thetaRadon)
        
        % Get the local maxima
        ind = locmax1d(averageRT(:,iTheta), winSize);
        ind = ind(averageRT(ind,iTheta) > 0);
        
        if numel(ind) ~= nnz(locMax(:,iTheta))
            ind
        end
        
        if ~isempty(ind)
            distAside = xp(ind);           
            distAlong = RD(ind,iTheta) ./ L(ind,iTheta);
            
            R = [ct(iTheta) st(iTheta); -st(iTheta) ct(iTheta)];
            pts = [distAside distAlong];
            
            % Store segment candidates parameters
            paramsCC = zeros(numel(distAside), 6);
            
            paramsCC(:,1:2) = repmat(cRadon - 1,numel(distAside),1) + pts * R;
            paramsCC(:,3) = averageRI(ind,iTheta);
            paramsCC(:,4) = sigmaPSF;
            paramsCC(:,5) = L(ind,iTheta);
            paramsCC(:,6) = theta(iTheta) + pi/2;

            % Generate the image model
            M = subResSegment2DImageModel(paramsCC,size(CC));
            
            % Assume segmentImageModel ~ N(mu,sigma^2) where mu is the mean
            % background intensity and sigma^2
            
            % Test the residual of each segment separatly (on the segment
            % support)
            % TODO
            
            % If resnorm is smaller than the one computed from a previous
            % angle, store the segments parameters
        end
    end
end
toc
%% --- Step 6 --- %%

%% --- Step 7 --- %%

%% --- Step 8 --- %%

params = [];
