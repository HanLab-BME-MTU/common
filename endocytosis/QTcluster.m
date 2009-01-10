function [clusterResults] = QTcluster(experiment,rest);

% QTcluster clusters initiations as defined by the rest vector in each
% movie in experiment and saves results to specified directory.
%
% INPUT:   experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the 
%                       function loadConditionData, and it needs to contain
%                       at least the field 
%                       .source, which is the (path) location of the 
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction                    
%           rest    =   restriction vector can have variable length;
%                       minimum length is five, where the entries are
%                       [stat da minfr minlft maxlft]
%                       optionally, the length can be extended to nine,
%                       where the additional entries are 
%                       [... minint maxint minmot maxmot]  
% OUTPUT
%           clusterResults.clusterResults = [xPos yPos clusterID lifetime startFrame]
%               the first column contains the x position, the second the y
%               position for all the initiations in a given movie
%               the third column contains an id number for each
%               particle, so that a particle at (x1,y1) belongs to cluster
%               that can be identified by the number in the third column of
%               row 1.
%           clusterResults.clusterCentroids = position of cluster centroids in two
%               columns; first column is x position and second column is y
%               position; 
%
% Uses: determineMovieLength
%       determineImagesize
%       RipleysKfunction
%       makeCorrFactorMatrix
%       QualityThresholdCluster
%
% Daniel Nunez, January 9, 2009

oldDir = cd;
%choose directory to save clustering results
[PATHNAME] = uigetdir(pwd, 'choose directory to save output');

%Fill in Missing Data
[experiment] = determineMovieLength(experiment);
[experiment] = determineImagesize(experiment);

%GET HOT SPOT RADIUS FROM DENSITY PLOTS
dist = 1:1:20;
for iexp = 1:length(experiment)
    
    %Load Lifetime Information
    cd([experiment(iexp).source filesep 'LifetimeInfo'])
    load('lftInfo')
    % status matrix
    statMat = lftInfo.Mat_status;
    % lifetime matrix
    lftMat = lftInfo.Mat_lifetime;
    % x-coordinate matrix
    matX = lftInfo.Mat_xcoord;
    % y-coordinate matrix
    matY = lftInfo.Mat_ycoord;
    % disapp status matrix
    daMat = lftInfo.Mat_disapp;
    % framerate
    framerate = experiment(iexp).framerate;
    % image size
    imsize  = experiment(iexp).imagesize;
    
    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
    
    msx = imsize(1);
    msy = imsize(2);
    imsizS = [imsize(2) imsize(1)];
    % construct convex hull out of complete point distribution
    % combined mpm
    selx = full(matX); selx = selx(isfinite(selx)); selx = nonzeros(selx(:));
    sely = full(matY); sely = sely(isfinite(sely)); sely = nonzeros(sely(:));
    combMPM = [ selx sely ];
    K = convhull(combMPM(:,1),combMPM(:,2));
    % edge points of the convex hull
    cpointsx = combMPM(K,1);
    cpointsy = combMPM(K,2);
    % create mask
    areamask = poly2mask(cpointsx,cpointsy,msx,msy);
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    corrFacMat = makeCorrFactorMatrix(imsizS, dist, 10, areamask'); 
    normArea = sum(areamask(:));

    mpmPits = [full(matX(findPos)) full(matY(findPos))];
    [kr,lr]=RipleysKfunction(mpmPits,mpmPits,imsizS,dist,[],normArea);
    Lpits(:,iexp) = lr;  
end

    [dlx,dly,dlz] = size(Lpits);
    carea = dist.^2;
    careadiff = carea; careadiff(2:length(careadiff)) = diff(carea);
    amat = repmat(careadiff',1,dly);
    dmat = repmat(dist',1,dly);
    currLR = Lpits;
    currKR = (currLR+dmat).^2;
    currKRdiff = currKR;
    currKRdiff(2:length(dist),:) = diff(currKR,1);
    currDen = currKRdiff./amat;
    pitDen = nanmedian(currDen,2);

%fit spline to productive data
ppProd = csaps(dist,pitDen);
denFitValuesProd = fnplt(ppProd);
%find max value of fit
maxNormDenProd = max(denFitValuesProd(2,:));
%find where radius is 10 pixels
findMidRadProd = find(denFitValuesProd(1,:) >= 10,1,'first');
%take the mean density from 10 pixels to the end as the far distance value
farDenProd = mean(denFitValuesProd(2,findMidRadProd:end));
%get point at 10 percent of difference in between maximum and long distance
%densities
findRadProd = find(denFitValuesProd(2,:) < 1+farDenProd,1,'first');
hotSpotRadiusProd = denFitValuesProd(1,findRadProd);
hotSpotRadius = hotSpotRadiusProd;


%%
%FOR EACH MOVIE
for iexp = 1:length(experiment)

    %Load Lifetime Information
    cd([experiment(iexp).source filesep 'LifetimeInfo'])
    load('lftInfo')
    % status matrix
    statMat = full(lftInfo.Mat_status);
    % lifetime matrix
    lftMat = full(lftInfo.Mat_lifetime);
    % x-coordinate matrix
    matX = full(lftInfo.Mat_xcoord);
    % y-coordinate matrix
    matY = full(lftInfo.Mat_ycoord);
    % disapp status matrix
    daMat = (lftInfo.Mat_disapp);
    % framerate
    framerate = experiment(iexp).framerate;
    % image size
    imsize  = experiment(iexp).imagesize;

    %CALCULATE NORMALIZED AREA OF CELL
    imagesize = experiment(iexp).imagesize;
    msx = imagesize(1);
    msy = imagesize(2);
    imsizS = [imagesize(2) imagesize(1)];
    % construct convex hull out of complete point distribution
    % combined mpm
    selx = full(matX); selx = selx(isfinite(selx)); selx = nonzeros(selx(:));
    sely = full(matY); sely = sely(isfinite(sely)); sely = nonzeros(sely(:));
    combMPM = [ selx sely ];
    K = convhull(combMPM(:,1),combMPM(:,2));
    % edge points of the convex hull
    cpointsx = combMPM(K,1);
    cpointsy = combMPM(K,2);
    % create mask
    areamask = poly2mask(cpointsx,cpointsy,msx,msy);
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    corfacmatM = makeCorrFactorMatrix(imsizS, dist, 10, areamask');
    normArea = sum(areamask(:));

    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
    alteredPos = findPos;
    %find pits that are too far away from other pits to be clustered
    alteredDistMat = squareform(pdist([matX(alteredPos) matY(alteredPos)]));
    %make zeros into nans since min of this distance matrix will always be
    %zero otherwise, which is the distance from one pit to its own self
    alteredDistMat(alteredDistMat == 0) = nan;
    findLonelyPits = find(min(alteredDistMat,[],2) > 2*hotSpotRadius);
    %store pits that are too far from other pits as pits outside hotspots
    outsidePits = alteredPos(findLonelyPits);
    %erase these pits from alteredPos
    alteredPos(findLonelyPits) = [];

    particlePositions = [matX(alteredPos) matY(alteredPos)];
    
    [clusteredParticles,clusterCentroids] = QualityThresholdCluster(particlePositions,2,hotSpotRadius);
    
    %add removed unclustered pits
    outsideParticles = [matX(outsidePits) matY(outsidePits) zeros(length(outsidePits),1)];
    clusteredParticles = [clusteredParticles; outsideParticles];
    
    %add lifetimes
    clusteredParticles(1:size(clusteredParticles,1),4) = lftMat([alteredPos;outsidePits])';
    %add start frame
    [dummy,startFrame] = ind2sub(size(lftMat),[alteredPos;outsidePits]);
    clusteredParticles(1:size(clusteredParticles,1),5) = startFrame';
    %make result structure for movie
    clusterResults(iexp).clusterResults = clusteredParticles;
    clusterResults(iexp).clusterCentroids = clusterCentroids;
    clusterResults(iexp).cellArea = normArea;
    clusterResults(iexp).hotSpotRadius = hotSpotRadius;
    clusterResults(iexp).framerate = framerate;
    clusterResults(iexp).movieID = iexp;
    clusterResults(iexp).movie = experiment(iexp).source;
end %for each movie
%save data onto folder under density folder
cd(PATHNAME);
mkdir(PATHNAME,'ClusterData')
filePath = [PATHNAME filesep 'ClusterData' filesep 'hotSpotAnalysisResultsWithMorePitsAnd1PlusProdFarDenRadius' datestr(now,'yyyymmdd')];
securesave(filePath,'clusterResults')
cd(oldDir)
end %of function
