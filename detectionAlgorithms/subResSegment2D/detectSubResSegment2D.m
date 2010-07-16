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
finalPoints = detectSubResFeatures2D(ima,cands,sigmaPSF,alpha,0,1,bitDepth,0,stdNoise);

%% --- Step 5 ---- %%

winSize = 2*ceil(2*sigmaPSF)+1;

thetaRadon = 0:4:179;
theta = -thetaRadon*pi/180;
ct = cos(theta);
st = sin(theta);
tt = tan(theta);
tic;

imagesc(segmentMask),colormap gray,axis image,axis off;

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
    
    % Crop R using bounding box and restrinct to CC's footprint
    filteredCC1 = zeros(bb(4),bb(3));
    filteredCC1(indLocal) = filteredIma1(CCstats(iCC).PixelIdxList);
    
    % Compute the Radon transform on the crop filtered image
    [RT xp] = radon(filteredCC1, thetaRadon);
    
    % Compute the Radon transform on the CC's footprint
    L = radon(CC, thetaRadon);
    
    % Compute the mean integral of R along the radon lines.
    averageRT = RT ./ (L + 1e-10);
    
    % The Radon origin is floor((size(BW)+1)/2).
    cRadon = floor((bb(3:4)+1)/2);
    
    hold on;
    plot(bb(1) + cRadon(1) - 1, bb(2) + cRadon(2) - 1,'rx');
    
    % Compute the Radon transform of a function defines as the signed
    % distance from the perpendical axis passing through cRadon. The
    % orientation of that perpendical follows each radon orientation.
    c1 = tt ./ sqrt(1 + tt.^2);
    c2 = 1 ./ sqrt(1 + tt.^2);
    D1 = zeros(size(CC));
    D1(indLocal) = cRadon(1) - x;
    D2 = zeros(size(CC));
    D2(indLocal) = y - cRadon(2);
    RD = bsxfun(@times, radon(D1, thetaRadon), c1) + bsxfun(@times, radon(D2, thetaRadon), c2);
    
    for iTheta = 40%1:numel(thetaRadon)
        
        % Get the local maxima
        ind = locmax1d(averageRT(:,iTheta), winSize);
        ind = ind(averageRT(ind,iTheta) > 0);
        
        if ~isempty(ind)
            distAside = xp(ind);           
            distAlong = RD(ind,iTheta) ./ L(ind,iTheta);
            
            R = [ct(iTheta) st(iTheta); -st(iTheta) ct(iTheta)];
            pts = [distAside distAlong];
            %pts = [distAside zeros(numel(ind), 1)];
            
            % Centers of Segments = center of the CC + R(pts)
            initialPos = repmat(bb(1:2) + cRadon - 1,numel(distAside),1) + pts * R;
                       
            % Store the length of each segment
            initialLength = L(ind,iTheta);
            
            % Store the orientation of each segment (note that theta is
            % perpendicular to the main segment's orientation. We add
            % pi/2 so the segment's orienation lies between [-pi/2 pi/2]
            initialAngle = repmat(theta(iTheta) + pi/2, numel(initialLength), 1);
            
            % TODO: check initialization
            
            % Compute optimization lsqnonlin
            % TODO
            
            % Test the residual of each segment separatly (on the segment
            % support)
            % TODO
            
            % If resnorm is smaller than the one computed from a previous
            % angle, store the segments parameters
            
            overlaySegment2DImage([], [initialPos zeros(size(initialLength)) initialLength initialAngle]);
        end
    end
end
toc
%% --- Step 6 --- %%

%% --- Step 7 --- %%

%% --- Step 8 --- %%
