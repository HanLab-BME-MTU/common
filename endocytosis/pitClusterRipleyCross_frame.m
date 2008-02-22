function [allResults] = pitClusterRipleyCross_frame(data, restrict, udist, scale);
% calculate the Ripley clustering and cross-clustering of two populations 
% defined by the restriction vector in each frame of the movie, then
% averaging over all frames and movies
% SYNOPSIS: [allResults] = pitClusterRipleyCross_frame(data, restrict, udist);
%
% INPUT:   data     = data structure, this should contain the fields
%                       .source
%                     which is the location where the TrackInfo is located,
%                     from which we can calculate the pertinent information
%                     like
%                       .lftInfo.Mat_status
%                       .lftInfo.Mat_xcoord
%                       .lftInfo.Mat_ycoord
%                       .lftInfo.Mat_disapp
%           restrict  = restriction vector can have variable length:
%           [stat da minfr minlft maxlft minint maxint minmot maxmot]
%           udist   =   distance vector for Ripley function
%           scale   =   distance scaling (optional), e.g. 0.067 (um/pix)
%
% OUTPUT:   allResults      = structure containing the fields:
%           .restriction    restrict vector for later reference
%           .dist          distance vector for later reference
%           .LR             all LR functions
%           .LRmean         LR averaged over all movies
%           .LRstd          std over all movies
%           .LRsem          standard error of the mean over all movies
%           .den            density of points of the two populations in
%                           each movie
%           
%
% written by Dinah Loerke
% last modified January 28th, 2008
% last modified Feb 13th, 2008

orDir = cd;

disp('extract MPMs from saved trackInfo');
% extract the relevant MPMs for population 1
[extractedMPMs1]=endoHotspotsMPMexDynamic(data, restrict(1,:));
% same for population 2
[extractedMPMs2]=endoHotspotsMPMexDynamic(data, restrict(2,:));

% loop over all available movies
for i=1:length(data)
    
    fprintf('movie #%02d',i);
    
    % relevant MPMs for this entry
    mpm1 = extractedMPMs1(i).mpm;
    [sx1,sy1] = size(mpm1);
    % number of frames
    nf1 = sy1/2;
    % number of defined points
    for k=1:nf1, devec = mpm1(:,2*k); np1(k) = length(find(devec>0)); end 
    
    mpm2 = extractedMPMs2(i).mpm;
    [sx2,sy2] = size(mpm2);
    nf2 = sy2/2;
    % number of points
    for k=1:nf2, devec = mpm2(:,2*k); np2(k) = length(find(devec>0)); end 
    
    % print the average number of points of each population in each frame -
    % if this number is extremely low, this acts as a cautionary in
    % interpreting the data
    print_np1 = round(nanmean(np1));
    print_np2 = round(nanmean(np2));
    fprintf(' nump1=%03d',print_np1);
    fprintf(' nump2=%03d',print_np2);
    
    % NOTE: the number of frames needs to be identical for both
    % extraction!!
    if nf1~=nf2
        error('number of frames doesn''t match!');
    else
        nf = nf1;
    end
        
        
    
    % image size
    imsiz = data(i).imagesize;
    msx = imsiz(1);
    msy = imsiz(2);
    imsizS = [imsiz(2) imsiz(1)];
    
       
    
    % determine current mask for correction
    % construct convex hull out of BOTH point distributions
    % combined mpm
    cx1 = mpm1(:,1:2:2*nf);
    cy1 = mpm1(:,2:2:2*nf);
    cx2 = mpm2(:,1:2:2*nf);
    cy2 = mpm2(:,2:2:2*nf);
    
    combMPMpoints = [ [cx1(:);cx2(:)] [cy1(:);cy2(:)] ];
    
    % defined positions
    defpos = find(combMPMpoints(:,1)>0);
    combMPMpoints = combMPMpoints(defpos,:);
    
    K = convhull(combMPMpoints(:,1),combMPMpoints(:,2));
    % edge points of the convex hull
    cpointsx = combMPMpoints(K,1);
    cpointsy = combMPMpoints(K,2);
    % IF DESIRED, create mask for correction here
    % note that the order is px,py,sy,sx
    areamask = [];
    areamask = poly2mask(cpointsx,cpointsy,msx,msy);  
         
        
   
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE
    % corfacmat = makeCorrFacMatrixGeneral(imsizeS, udist, 10); 
    corfacmatM = makeCorrFactorMatrix(imsizS, udist, 10, areamask'); 
    normArea = sum(areamask(:));
    
                    
    % self-correlation pop1
    [kr,lr11]=RipleysKfunction(mpm1, mpm1, imsizS, udist, corfacmatM,normArea);

    % self-correlation pop2
    [kr,lr22]=RipleysKfunction(mpm2, mpm2, imsizS, udist, corfacmatM,normArea);

    % cross-correlation pop1-pop2
    [kr,lr12]=RipleysKfunction(mpm1, mpm2, imsizS, udist, corfacmatM,normArea);
    % cross-correlation pop2-pop1
    % should be identical except for edge effects
    %[kr,lr21]=RipleysKfunction(mpm2, mpm1, imsizS, udist, corfacmatM,normArea);
                      
    lrAll(:,1) = nanmean(lr11,2);
    lrAll(:,2) = nanmean(lr22,2);
    lrAll(:,3) = nanmean(lr12,2);
    %lrAll(:,4) = nanmean(lr21,2);
    
    allLR(:,:,i) = squeeze(lrAll);
    allDen(1,i) = 1000 * nanmean(np1)/normArea;
    allDen(2,i) = 1000 * nanmean(np2)/normArea;
    
    fprintf('\n');
    
end % of for-loop

lrAV = nanmean(allLR,3);
lrSTD = nanstd(allLR,1,3);
lrSEM = lrSTD/sqrt(length(data));

allResults.restriction = restrict;
allResults.dist = udist;
allResults.LR = allLR;
allResults.LRmean = lrAV;
allResults.LRstd = lrSTD;
allResults.LRsem = lrSEM;
allResults.dens = allDen;
    
    

areaplot1 = [ ((lrAV(:,1)-lrSEM(:,1))') ; 2*lrSEM(:,1)' ];
areaplot2 = [ ((lrAV(:,2)-lrSEM(:,2))') ; 2*lrSEM(:,2)' ];


pdist = udist;
if nargin>3
    pdist = scale*udist;
end

figure
hold on;
area(pdist,areaplot1');
area(pdist,areaplot2');
plot(pdist,lrAV(:,1),'g-');
plot(pdist,lrAV(:,2),'r-');


%%=========================================================================
%           convert allLR to local densities
%==========================================================================


allDen = allLR;

carea = udist.^2;
careadiff = carea; careadiff(2:length(careadiff)) = diff(carea);
amat = repmat(careadiff',1,3);

dmat = repmat(udist',1,3);

for i=1:length(data)
    currLR = allLR(:,:,i);
    currKR = (currLR+dmat).^2;
    currKRdiff = currKR;
    currKRdiff(2:length(udist),:) = diff(currKR,1);
    currDen = currKRdiff./amat;
    
    allDen(:,:,i) = currDen;
end


ndenAV = nanmean(allDen,3);
ndenSTD = nanstd(allDen,1,3);
ndenSEM = ndenSTD/sqrt(length(data));

allResults.ND = allDen;
allResults.NDmean = ndenAV;
allResults.NDstd = ndenSTD;
allResults.NDsem = ndenSEM;



areaplot1d = [ ((ndenAV(:,1)-ndenSEM(:,1))') ; 2*ndenSEM(:,1)' ];
areaplot2d = [ ((ndenAV(:,2)-ndenSEM(:,2))') ; 2*ndenSEM(:,2)' ];

figure
hold on;
area(pdist,areaplot1d');
area(pdist,areaplot2d');
plot(pdist,ndenAV(:,1),'g-');
plot(pdist,ndenAV(:,2),'r-');


end % of function