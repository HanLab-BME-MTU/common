function [allResults]=pitClusterRipleyCross_project(data, restrict, udist, scale)
% calculate the Ripley clustering and cross-clustering of two populations 
% defined by the restriction vector projected over all frames of the
% movie, then averaging over all movies
% SYNOPSIS: [allResults]=pitClusterRipleyCross_project(data, udist, restrict)
%
% INPUT:   data         =   data structure, this should contain the fields
%                           .source
%                           which is the location where the TrackInfo is 
%                           located, from which we can calculate the 
%                           pertinent information like
%                           .lftInfo.Mat_status
%                           .lftInfo.Mat_xcoord
%                           .lftInfo.Mat_ycoord
%                           .lftInfo.Mat_disapp
%           restrict    =   restriction vector can have variable length:
%                           [stat da minfr minlft maxlft ...
%                               minint maxint minmot maxmot]
%           udist       =   distance vector for Ripley function
%           scale (opt) =   pixel sclaing, e.g. 0.067 if 1 pixel = 67nm
%
% OUTPUT:   allResults      = structure containing the fields:
%           .restriction    restrict vector for later reference
%           .dist          distance vector for later reference
%           .LR             all LR functions
%           .LRmean         LR averaged over all movies
%           .LRstd          std over all movies
%           .LRsem          standard error of the mean over all movies
%           .dens            density of points of the two populations in
%                           each movie
%
% written by Dinah Loerke
% last modified January 28th, 2008
% last modified Feb 13th, 2008


% number of restriction 'sub-populations'
[ccx,ccy] = size(restrict);
% total number of self- and cross-correlations
numcc = sum([1:ccx]);

% allLR
allLR = zeros(length(udist),numcc,length(data));


% loop over all fields in the data
for i=1:length(data)
    
    fprintf('movie #%02d',i);
    
    % current lftInfo 
    currLI = data(i).lftInfo;
    % status matrix
    mat_stat = currLI.Mat_status;
    % lifetime matrix
    mat_lft = currLI.Mat_lifetime;
    % x-coordinate matrix
    mat_x = currLI.Mat_xcoord;
    % y-coordinate matrix
    mat_y = currLI.Mat_ycoord;
    % disapp status matrix
    mat_da = currLI.Mat_disapp;
    % framerate
    fr = data(i).framerate;
    % image size
    imsiz = data(i).imagesize;
    msx = imsiz(1);
    msy = imsiz(2);
    imsizS = [imsiz(2) imsiz(1)];
       
    
    % construct convex hull out of complete point distribution
    % combined mpm
    selx = full(mat_x); selx = selx(isfinite(selx)); selx = nonzeros(selx(:));
    sely = full(mat_y); sely = sely(isfinite(sely)); sely = nonzeros(sely(:));
    combMPM = [ selx sely ];
    K = convhull(combMPM(:,1),combMPM(:,2));
    % edge points of the convex hull
    cpointsx = combMPM(K,1);
    cpointsy = combMPM(K,2);
    % create mask
    areamask = poly2mask(cpointsx,cpointsy,msx,msy);
    
    
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    corfacmatM = makeCorrFactorMatrix(imsizS, udist, 10, areamask'); 
    normArea = sum(areamask(:));
    
    cc_counter = 1;
    
    % loop over all rows of the restriction values, which designate
    % different conditions
    
    for r=1:ccx
        
        % for the current value of restvector, collect those locations of 
        % points that fulfill the specified restrictions, regarding
        % appearance/disappearance status and minimum/maximum lifetime

        %desired minimum lifetime in seconds
        dstat   = restrict(r,1);
        dapp    = restrict(r,2);
        dminfr  = restrict(r,3);
        minlft  = restrict(r,4);
        maxlft  = restrict(r,5);
        minlft_fr = round(minlft/fr);
        maxlft_fr = round(maxlft/fr);

        [kx,ky]=size(mat_da);
    
        findpos = find( (mat_stat==dstat) & (mat_da==dapp) &...
            (mat_lft>dminfr) & (mat_lft>minlft_fr) & (mat_lft<maxlft_fr));
        findx = full(mat_x(findpos));
        findy = full(mat_y(findpos));
        
        currMPM = [findx findy];
        restrictedMPMs(r).mpm = currMPM;
        
        % density of objects
        allDen(r,i) = 1000 * length(findpos)/sum(areamask(:));
        
         % print the average number of points of each population in each 
         % frame - if this number is extremely low, this acts as a 
         % cautionary in interpreting the data
        print_np = length(findpos);
        fprintf([' nump',num2str(r)]);
        fprintf('=%05d',print_np);
            
    end
    
   
    
    
    % now loop over all combinations of the different restriction
    % conditions
    for k1=1:ccx
        
        currMPM1 = restrictedMPMs(k1).mpm;
                        
        % determine second mpm file
        for k2 = k1:ccx
            
            currMPM2 = restrictedMPMs(k2).mpm;
            
                       
            % determine mask for correction
            if (length(currMPM1)>20) & (length(currMPM2)>20)
                              
                % determine ripley clustering function for
                % cross-correlation of the two specified mpms (which may be
                % the same mpm in some cases)       

                % cross-correlation mpm1-mpm2
                %[krCross,lrCross]=RipleyKfunc_standalone(currMPM1, currMPM2, imsizS, udist, corfacmat);
                [krCross,lrCross]=RipleysKfunction(currMPM1, currMPM2, imsizS, udist, corfacmatM, normArea);
                
                % store for averaging
                allLR(:,cc_counter,i) = lrCross;

                                
                % cross_correlation combination
                ccComb(cc_counter,:) = [k1,k2];

                
                % display
                % hold on; plot(udist,lrCross,'b.-'); 
                % hold on; plot(udist,lrCross2,'r.-');
                % pause(0.1);
                
            else
                allLR(:,cc_counter,i) = nan;
                allDen(i) = nan;
            end % of if
            
            cc_counter = cc_counter+1;
            
                      
            
        end % of for k2
        
    end % of for k1
            
    fprintf('\n');
    
end




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
allResults.comb = ccComb;

 
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
       