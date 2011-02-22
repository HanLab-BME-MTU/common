function [clusters, ImBG] = detectSegment2D(ima,mask,sigmaPSF,minSize,bitDepth,display)

if nargin <= 5 || isempty(display)
    display = false;
end

%% --- Initialization --- %%

% Make sure image's class is double
if ~isa(ima,'double')
    ima = double(ima);
end

[nrows,ncols] = size(ima);

% Parameters for detectSubResFeature2D

% %alpha.alphaR = .05;
% alpha.alphaA = .01;
% % alpha.alphaD = .1;
% alpha.alphaF = 0;

% Parameters for segment2DFit

options = optimset('Jacobian', 'on', ...
    'MaxFunEvals', 500, ...
    'MaxIter', 500, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

%% --- Step 1: Create Segment Mask ---- %%

% Filter image with 2nd order steerable filter
% TODO: Use M=4 instead
filteredIma1 = steerableFiltering(ima,2,sigmaPSF);
filteredIma1(filteredIma1 < 0) = 0;
filteredIma1(mask == false) = 0;

% Segment filteredIma1
segmentMask = logical(blobSegmentThreshold(filteredIma1,minSize,0,mask));

%% --- Step 2: Create Point Mask ---- %%

% Filter image with 2-D laplacian filter
filteredIma2 = filterLoG(ima,sigmaPSF);
filteredIma2(~mask) = 0;

% Filter image with a Gaussian filter to compute bkg value.
filteredIma3 = Gauss2D(ima,sigmaPSF);
filteredIma3(~mask) = 0;

% Get local min/max of filteredIma2.
hside = 3 * ceil(sigmaPSF) + 1;
locMax = locmax2d(filteredIma2, [hside hside]);
locMin = locmin2d(filteredIma2, [hside hside]);

indMin = find(locMin);
indMax = find(locMax);
[yMin xMin] = ind2sub(size(locMin), indMin);
[yMax xMax] = ind2sub(size(locMax), indMax);
tri = delaunay(xMin,yMin);
triangles = tsearch(xMin,yMin,tri,xMax,yMax);
nnzIdx = ~isnan(triangles);
tri2 = tri(triangles(nnzIdx),:);

avgBkg = inf(size(indMax));
avgBkg(nnzIdx) = mean(reshape(filteredIma3(indMin(tri2(:))), size(tri2)),2);
stdNoise = std(ima(indMin));

pointCands = zeros(size(locMax));
indMax = indMax(ima(indMax) > avgBkg + 2 * stdNoise);
pointCands(indMax) = filteredIma2(indMax);
 
% figure
% imshow(ima,[])
% hold on;
% ind = find(locMax);
% [y x] = ind2sub(size(pointCands),ind);
% plot(x,y,'b.')
% ind = find(pointCands);
% [y x] = ind2sub(size(pointCands),ind);
% plot(x,y,'g.')
% nnz(indMax)

% DEPRECATED
% Generate pointMask from pointCands
% pointDist = bwdist(pointCands ~= 0);
% pointMask = pointDist < hside;
% 
% % Make the union of the segment and point masks
% unionMask = segmentMask | pointMask;

% Compute property of connected components
CCstats = regionprops(segmentMask,'PixelIdxList','BoundingBox');
nCC = numel(CCstats);

% Define the resulting set of clusters
clusters(1:nCC) = struct('avgBkg',[],'stdBkg',[],'segments',[],'points',[]);

% Create background image (original image without feature areas)
bkgImage = ima;
bkgImage(segmentMask) = 0;

%% --- Step 4 ---- %%

% TODO: Possible optimization: instead of computed Radon transform of
% filteredIma1 for every possible angle, use the Feature-adapted Radon
% transform.

winSize = 2*ceil(2*sigmaPSF)+1;

thetaRadon = 0:179;
theta = -thetaRadon*pi/180;
ct = cos(theta);
st = sin(theta);
tic;
for iCC = 1:nCC

    % Floor bounding box
    bb = ceil(CCstats(iCC).BoundingBox);
    
    % Get the footprint of the CC
    CC = zeros(bb(4),bb(3));
    [y x] = ind2sub([nrows ncols],CCstats(iCC).PixelIdxList);
    x = x - bb(1) + 1;
    y = y - bb(2) + 1;
    indLocal = sub2ind(size(CC), y, x);
    areaCC = numel(indLocal);
    CC(indLocal) = 1;
    
    %% Compute initial background value of the cluster
    
    % we use a 1-pixel extended bounding box to have sufficient background
    % pixels.
    bkgImageCC = bkgImage(max(bb(2)-1,1):min(bb(2)+bb(4),nrows), ...
        max(bb(1)-1,1):min(bb(1)+bb(3),ncols));
    
    bkgPixels = nonzeros(bkgImageCC);
    
    clusters(iCC).avgBkg = mean(bkgPixels);
    clusters(iCC).stdBkg = std(bkgPixels);
    
    %%
    
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

    % TODO: We can even make more tests here to discard unsignificant local
    % maxima.
    
    %% Compute initial value of the center's coordinates of each segment candidates
    
    % Compute the distance of each local maxima away from the main axis,
    % for each main axis's orientation
    distAside = bsxfun(@times,locMax,xp);
    
    % Compute the distance of each local maxima away form the secondary
    % axis, for each secondary axis's orientation. For that, compute the
    % Radon transform of a function defines as the signed distance from the
    % perpendical axis passing through cRadon. The orientation of that
    % perpendical follows each radon orientation.plot(pointCands
    D1 = zeros(size(CC));
    D1(indLocal) = cRadon(1) - x;
    D2 = zeros(size(CC));
    D2(indLocal) = y - cRadon(2);
    RD = bsxfun(@times, radon(D1, thetaRadon), st) + ...
        bsxfun(@times, radon(D2, thetaRadon), ct);
    
    distAlong = zeros(size(locMax));
    distAlong(locMax) = RD(locMax) ./ L(locMax);
    
    xCoord = bsxfun(@times,distAside,ct) + bsxfun(@times,distAlong,-st);
    xCoord(locMax) = cRadon(1) + xCoord(locMax);
    yCoord = bsxfun(@times,distAside,st) + bsxfun(@times,distAlong, ct);
    yCoord(locMax) = cRadon(2) + yCoord(locMax);
    
    %% Compute initial amplitude of each segment
    
    % Crop ima
    imaCC = ima(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1);
    
    % Compute the Radon transform on the crop image. averageRI correspond
    % to the mean intensity of the image along the radon lines.
    RI = radon(imaCC, thetaRadon);
    
    amp = zeros(size(locMax));
    amp(locMax) = RI(locMax) ./ L(locMax);

    %% Compute initial value of the length of each segment candidates
    
    length = zeros(size(locMax));
    length(locMax) = L(locMax);
    
    %% Point classification
    
    % To choose the right segment's orientation, we classify the
    % local maxima detected in STEP 3 and choose the orientation that
    % yields the largest likelihood. We follow the paper "Finding
    % Curvilinear Features in Spatial Point Patterns: Principal Curve
    % Clustering with Noise", IEEE Transactions on PAMI, Vol 22, No 6,
    % 2000. Derek C. and Adrian E. Raftery.
    
    % Crop the image of local maxima to the CC    
    pointCandsCC = zeros(bb(4),bb(3));
    pointCandsCC(indLocal) = pointCands(CCstats(iCC).PixelIdxList);
    
    % Find the point candidates
    indPointCandsCC = find(pointCandsCC);
    nPointCandsCC = numel(indPointCandsCC);
    [y x] = ind2sub(size(pointCandsCC), indPointCandsCC);   
    
    % Store BIC for each model (i.e. each orientation)
    bic = Inf(numel(theta),1);
    
    pointClasses = zeros(numel(theta), nPointCandsCC);
    
    % Empty class
    c0 = repmat(1/areaCC,1, nPointCandsCC);
    
    % For each orientation...

    for iTheta = 1:numel(theta)
        
        ind = find(xCoord(:,iTheta));
        
        pX = bsxfun(@minus,xCoord(ind,iTheta),x');
        pY = bsxfun(@minus,yCoord(ind,iTheta),y');

        distAbout = bsxfun(@times,ct(iTheta), pX) + bsxfun(@times,st(iTheta), pY);
        
        likelihood = 1 / (sqrt(2 * pi) * sigmaPSF) * exp(-.5 * distAbout.^2 / sigmaPSF);
        
        [~, I] = max([likelihood; c0],[],1);
        
        pointClasses(iTheta,I ~= numel(ind) + 1) = ind(I(I ~= numel(ind) + 1));
        
        % Number of model parameters
        k = 3 * numel(ind) + 1;
        
        bic(iTheta) = - 2 * log(prod(sum(likelihood,1))) + k * log(nPointCandsCC);
    end
    
    % Find the most likely set of segment candidates
    [~, iTheta] = min(bic);
    ind = find(locMax(:,iTheta));
    nSegments = numel(ind);
    
    %% --- STEP 5: Optimize segment candidates --- %%
    
    if nSegments

        clusters(iCC).segments = horzcat(xCoord(ind,iTheta), ...
            yCoord(ind,iTheta), ...
            amp(ind,iTheta), ...
            repmat(sigmaPSF, nSegments, 1), ...
            length(ind,iTheta), ...
            repmat(theta(iTheta) + pi/2, nSegments, 1));
    
        % optimize the length and amplitude
        paramSelector = false(1,6);
        paramSelector([3 5]) = true;
        
        varSegmentParams = clusters(iCC).segments(:,paramSelector);
        fixSegmentParams = clusters(iCC).segments(:,~paramSelector);
        
        fun = @(params) segment2DFit(params, imaCC, clusters(iCC).avgBkg, ...
            fixSegmentParams, paramSelector);
        
        lb = [zeros(nSegments,1) varSegmentParams(:,1) - 3 * clusters(iCC).stdBkg];
        ub = varSegmentParams;
        ub(:,2) = ub(:,2) + 3 * clusters(iCC).stdBkg;
        
        %disp(num2str(iCC));
        
        varSegmentParams = lsqnonlin(fun, varSegmentParams,ub,lb,options);
        
        clusters(iCC).segments(:,paramSelector) = varSegmentParams;
        
        clusters(iCC).segments(:,1:2) = clusters(iCC).segments(:,1:2) + ...
            repmat(bb(1:2),numel(ind),1) - 1;    

        % Remove segment whose length <= 1
        %clusters(iCC).segments = clusters(iCC).segments(clusters(iCC).segments(:,5) > minSize,:);
    end
    
    %% --- STEP 6: Optimize point candidates --- %%
    
    % Find the point candidates
    ind = find(pointClasses(iTheta,:) == 0);
    
    if ~isempty(ind)
        pts = [x(ind) y(ind)];
%         
%         % Convert points to speckle cands to comply with Khuloud format.
%         nPSF = numel(ind);
%         cands(1:nPSF) = struct('Lmax',[],'IBkg',[],'status',[]);
%         
%         [cands(:).IBkg] = deal(clusters(iCC).avgBkg);
%         [cands(:).status] = deal(true);
%         
%         for iPSF = 1:nPSF
%             cands(iPSF).Lmax = pts(iPSF,[2 1]);
%         end
%                 
% %         subResPts =
% detectSubResFeatures2D(imaCC,cands,sigmaPSF,alpha,0,1, ...
% %             bitDepth,0,clusters(iCC).stdBkg);
% %         
% %         % Store the final set of independent points.
% %         clusters(iCC).points = [subResPts.xCoord' subResPts.yCoord'] + ...
% %             repmat(bb(1:2),numel(subResPts.xCoord),1) - 1;
clusters(iCC).points = pts + repmat(bb(1:2),size(pts,1),1) - 1;
     end
end
toc

%% DISPLAY
if display
    S = cell2mat(arrayfun(@(iCC) clusters(iCC).segments, (1:nCC)', 'UniformOutput', false));
    P = cell2mat(arrayfun(@(iCC) clusters(iCC).points, (1:nCC)', 'UniformOutput', false));
    
    if ~isempty(S)
        overlaySegment2DImage(ima,S);
    else
        imshow(ima,[]);
    end
    
    hold on;
    
    if ~isempty(P)
        plot(P(:,1),P(:,2),'r.');
    end
    
    % Display point
    ind = find(pointCands);
    [y x] = ind2sub(size(pointCands),ind);
    plot(x,y,'r.');
    
    hold off;
end
