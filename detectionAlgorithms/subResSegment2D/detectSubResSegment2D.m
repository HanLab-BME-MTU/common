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

finalPoints = detectSubResFeatures2D(ima,cands,sigmaPSF,alpha,0,1,bitDepth,0,stdNoise);

%% --- Step 5 ---- %%

thetaRadon = 0:179;
theta = -thetaRadon*pi/180;
ct = cos(theta);
st = sin(theta);
tt = tan(theta);

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
    
    % Compute the radon transform on the crop filtered image
    [RT xp] = radon(filteredCC1, thetaRadon);
    
    % Compute the radon transform on the CC's footprint
    L = radon(CC);
    
    % Compute the mean integral of R along the line oriented along t90.
    averageRT = RT ./ (L + 1e-10);
    
    % The Radon origin is floor((size(BW)+1)/2).
    cRadon = floor((bb(:,3:4)+1)/2);
    
    for iTheta = 1:numel(thetaRadon)
        
        % Get the local maxima
        ind = locmax1d(averageRT(:,iTheta), winSize);
        ind = ind(averageRT(ind) > 0);
        
        if ~isempty(ind)
            distAside = xp(ind);
            
            % D is the signed distance transform from the perpendicular axis.
            D = zeros(bb(4),bb(3));
            D(indLocal) = (tt(iTheta) * (cRadon(1) - x) + y - cRadon(2)) ./ ...
                sqrt(1 + tt(iTheta)^2);
    
            % The integration of that distance transform along lines and restricted
            % to the CC's footprint will gives the center of each FA.
            distAlong = radon(D,thetaRadon);
            distAlong = distAlong(indMax) ./ L;
            
            R = [ct(iTheta) st(iTheta); -st(iTheta) ct(iTheta)];
            pts = [distAside distAlong];
            
            % Centers of Segments = center of the CC + Rot90(pts)
            initialPos = repmat(bb(:,1:2) + cRadon - 1,numel(distAside),1) + pts * R;
                       
            % Store the length of each segment
            initialLength = L(ind);
            
            % Store the orientation of each segment (note that theta is
            % perpendicular to the main segment's orientation. We add
            % pi/2 so the segment's orienation lies between [-pi/2 pi/2]
            initialAngle = repmat(theta(iTheta) + pi/2, numel(initiPos), 1);
            
            % TODO: check initialization
            
            % Compute optimization lsqnonlin
            % TODO
            
            % Test the residual of each segment separatly (on the segment
            % support)
            % TODO
            
            % If resnorm is smaller than the one computed from a previous
            % angle, store the segments parameters
        end
    end
end

%% --- Step 6 --- %%

%% --- Step 7 --- %%

%% --- Step 8 --- %%
