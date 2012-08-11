function [kr,lr,pcr,glr]=RipleysKfunction(mpm1,mpm2,imsiz,dist,corrFacMat,normArea,sparse)
% RipleysKfunction calculates Ripley's K-function for a given MPM,
% allowing cross-corrlation between two MPMs
% SYNOPSIS  [kr,lr,pcr]=RipleysKfunction(mpm,imsiz,dist,corrFacMat, normArea);
%
% INPUT     mpm1:      mpm file containing (x,y) coordinates of points in
%                      the image in succesive columns for different time
%                      points
%           mpm2:      second mpm
%           imsiz:     x,y-size of the image (maximum possible value for x-coordinate)
%           dist:      distance vector e.g. [1:50]
%           corrMat:   OPTIONAL if a correction matrix is pre-calculated
%                      outside of this function, it can be used directly
%           normArea:  OPTIONAL if the point density for normalization
%                      should not be based on number of points per total
%                      rectangular image size, but e.g. per a different
%                      area of interest (a mask of which may have already
%                      been used to calculate corrMat), then the area
%                      for normalization should be entered here
%           sparse: true to measure distances using graph based algorithm
%                   that creates sparse matrix (might be faster for large
%                   data sets)
%
%           NOTE: IF you want to pre-calculate the correction factor matrix
%           (which is recommended because it saves time), then run the
%           function:
%           [corrFacMatrix] = makeCorrFactorMatrix(imsiz, dist, samplesize,
%           mask);
%
%           IF nargin == 5              the provided corrFacMat is used
%           IF corrFacMat is an integer corrFacMat is calculated using
%                                       this integer value as samplesize
%           IF nargin == 4              Ripley's edge correction is used
%           IF nargin == 3 (or ==nan)   dist is set as default [1:rs]
%
% OUTPUT
%           kr:     Ripley's K-function
%           lr:     Besag's L-function = sqrt(K(r))-r
%           pcr:    pair-correlation function
%           glr:     Getis' L-function = sqrt(K(r))
%           for all, every column contains the kr/lr function for one
%           frame of the mpm-file, the row values correspond to the
%           specified distances
%
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
%
% last modified
% DATE:     29-Jan-2008 (last update)
%           1-July-2009
%
%

% create vector containing x- and y-image size
imsizex = imsiz(1);
imsizey = imsiz(2);

% if no distance vector is chosen explicitly or if dist is nan, the
% distance vector is set as default to 1:rs, where rs is the half diagonal
% of the image (this is the standard in the literature)
rs = round(sqrt(imsizex^2+imsizey^2));
if nargin<4
    distvec = 1:rs;
    nr = rs;
elseif isnan(dist)
    distvec = 1:rs;
    nr = rs;
else
    distvec = dist;
    nr = length(distvec);
end

% if corrFacMat is specified but consists only of an integer value,
% calculate the matrix using that value as sample, else set corrFacMat to
% empty
if (nargin<5)
    corrFacMat = [];
else
    if length(corrFacMat)==1
        sample = corrFacMat;
        corrFacMat = makeCorrFactorMatrix(imsiz,sample);
    end
end


if nargin<7
    sparse = false;
end

%determine size of mpm-file
[nx1,ny1]=size(mpm1);
%number of frames
numframes1 = round(ny1/2);
%determine size of mpm-file
[nx2,ny2]=size(mpm2);
%number of frames
numframes2 = round(ny2/2);

% function requires both mpms to have the same number of frames
if numframes1~=numframes2
    error('number of frames in 2 mpms doesn''t match');
else
    numf = numframes1;
end

%initialize results matrix pvr; x-dimension equals the employed number of
%values for the circle radius, y-dimension equals number of planes of the
%input mpm-file
kr  = zeros(nr,numf);
lr  = kr;
pcr = lr;
glr  = lr;

% loop over all frames
for i=1:numf
    % relevant nonzero x,y-positions in mpm1
    pos1    = find(mpm1(:,i*2)>0);
    cmpm1   = mpm1(pos1,2*i-1:2*i);
    np1     = length(pos1);
    
    % relevant nonzero x,y-positions in mpm2
    pos2    = find(mpm2(:,i*2)>0);
    cmpm2   = mpm2(pos2,2*i-1:2*i);
    np2     = length(pos2);
    
    %fprintf(' frame %04d',i);
    
    % if there are any relevant points in this frame - at least one point
    % each for cross-correlation (different matrices mpm1/mpm2), at least
    % two points for auto-corr (mpm1 == mpm2)
    if ( min(np1,np2)>0 ) && ( (np1+np2)>2 )
        
        % output pr: # of points as a function of distance
        [pr,nump] = pointsincircleCross(cmpm1,cmpm2,imsiz,distvec,corrFacMat,sparse);
        
        % normalized pr - normalize by total point density
        totaldensity = (nump)/(imsizex*imsizey);
        if nargin>5 && ~isempty(normArea)
            totaldensity = nump/normArea;
        end
        
        prnorm = (pr/pi)*(1/totaldensity);
        %prnorm = pr;
        kr(:,i)     = prnorm;
        lr(:,i)     = sqrt(prnorm) - distvec;
        pcr(:,i)    = convertLR2PCF(lr(:,i),distvec);
        glr(:,i) = sqrt(prnorm);
    else
        kr(:,i)     = nan*dist;
        lr(:,i)     = nan*dist;
        pcr(:,i)    = nan*dist;
        glr(:,i)    = nan*dist;
    end
    
    %fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for i-loop

%fprintf('\n');

end % of function





%=========================================================================
%=========================================================================
%=========================================================================
%====================       SUBFUNCTIONS    ==============================
%=========================================================================
%=========================================================================
%=========================================================================




function [npvr,nump]=pointsincircleCross(m1,m2,ms,dist,corrMat,sparse)
% pointsincircle calculates the average number of points in a circle around
% a given point as a function of the circle radius (averaged over all points
% and normalized by total point density); this function is called Ripley's
% K-function in statistics, and is an indication of the amount of clustering
% in the point distribution
%
% SYNOPSIS   [m2,num]=pointsincircle(m1,m2,ms,dist,corrMat);
%
% INPUT     m1:     matrix of size (n1 x 2) containing the (x,y)-coordinates
%                   of n1 points; these points are considered the CHILDREN
%           m2:     matrix of size (n2 x 2) containing the (x,y)-coordinates
%                   of n2 points; these points are considered the PARENTS
%           ms:     vector containing the parameters [imsizex imsizey] (the
%                   x-size and y-size of the image)
%           dist:   distance vector
%           corrMat:   correctionFactor matrix; if this matrix is empty or
%                   nargin<5, then the simple Ripley's edge correction is
%                   used to calculate the correction factor
%           sparse: true to measure distances using graph based algorithm
%                   that creates sparse matrix (might be faster for large
%                   data sets)
%
%
% OUTPUT    npvr:   vector containing the number of points in a circle
%                   around each point, for an increasing radius;
%                   radius default values are 1,2,3,....,min(ms)
%                   function is averaged over all objects in the image
%           nump:   number of points
%
% Dinah Loerke, Jan 29th, 2008

nump1 = size(m1,1);
msx=ms(1);
msy=ms(2);

% in the following, the matrix mdist will contain the distance of all
% points in m1 from all points in m2, where m1 are the children and m2 are
% the parents; because of the way the function below is set up in terms of
% rows-columns, the order in the DistanceMatrix function needs to be
% parent-child

%Distance Matrix is faster for smaller matrices, but
%createSparseDistanceMatrix can better deal with large matrices
if sparse
    mdist = createSparseDistanceMatrix([m2(:,2) m2(:,1)],[m1(:,2) m1(:,1)],max(dist(:)),0);
else
    [mdist]=DistanceMatrix(m2,m1);
end

% NOTE: In the old version of the Ripley, since it was designed for
% self-correlation, only distances > 0 were considered in the distance
% histogram below. In this version, since it allows cross-correlation,
% zero-distances have to be included for DIFFERENT mpms, but should be
% excluded for IDENTICAL mpms

% monitor progress
%fprintf(' progress %02d',0);

if isempty(m1) || isempty(m2) || isempty(mdist)
    npvr = nan*dist;
    nump = 0;
else
    
    %allocate space
    histmat = nan(length(mdist(:,1)),length(dist)+1);
    corrmat = nan(length(mdist(:,1)),length(dist));
    for i=1:length(mdist(:,1))
        
        % matrix histmat contains the distance histogram for each point, where
        % every row represents the distance histogram for the cell at that
        % position
        histmat(i,:) = histc(nonzeros(mdist(i,:)),[0 dist]);
        
        % determine the circumference correction vector for this point; in the
        % matrix corrmat, the columns represent the points, and the rows
        % represent the entries for the radius vector dist
        % IF corrMat exists, use these values directly, else calculate corrmat
        % in situ with Ripley's edge correction
        if ~isempty(corrMat)
            cpx = round(m2(i,1)); cpy = round(m2(i,2));
            corrmat(i,:)= corrMat(max(cpx,1),max(cpy,1),:);
        else
            % x- and y-distances of this center point (the parent point) from the
            % nearest edge
            ex = min(m2(i,1),1+msx-m2(i,1));
            ey = min(m2(i,2),1+msy-m2(i,2));
            corrmat(i,:)= circumferenceCorrFactor(ex,ey,dist,msx,msy);
        end
        
        %         % update iter
        %         iter = round( 100*(i/length(mdist(:,1))) );
        %         %if iter<100, fprintf('\b\b%02d',iter); end
        
    end
    
    
    
    % multiply the matrices histmat and the correction vector
    histmat(:,length(dist)+1)=[];
    pointsmat = histmat./corrmat;
    
    % for every point, the (number-of-points within radius r)-function is the
    % cumulative sum of the values in the corresponding row
    pointsincirclemat = cumsum(pointsmat,2);
    
    % average for this frame over all existing points
    npvr = nanmean(pointsincirclemat,1);
    nump = nump1;
    
end % of if there are any points



end % of function


function [m2]=DistanceMatrix(c1,c2)
%this subfunction makes a neighbour-distance matrix for input matrix c1
%(n1 x 2 points) and c2
%output: m2 (n1 x n1) matrix containing the distances of each point in c1
%from each point in c2

np1=size(c1,1);function [kr,lr,pcr,glr]=RipleysKfunction(mpm1,mpm2,imsiz,dist,corrFacMat,normArea,sparse)
% RipleysKfunction calculates Ripley's K-function for a given MPM,
% allowing cross-corrlation between two MPMs
% SYNOPSIS  [kr,lr,pcr]=RipleysKfunction(mpm,imsiz,dist,corrFacMat, normArea);
%
% INPUT     mpm1:      mpm file containing (x,y) coordinates of points in
%                      the image in succesive columns for different time
%                      points
%           mpm2:      second mpm
%           imsiz:     x,y-size of the image (maximum possible value for x-coordinate)
%           dist:      distance vector e.g. [1:50]
%           corrMat:   OPTIONAL if a correction matrix is pre-calculated
%                      outside of this function, it can be used directly
%           normArea:  OPTIONAL if the point density for normalization
%                      should not be based on number of points per total
%                      rectangular image size, but e.g. per a different
%                      area of interest (a mask of which may have already
%                      been used to calculate corrMat), then the area
%                      for normalization should be entered here
%           sparse: true to measure distances using graph based algorithm
%                   that creates sparse matrix (might be faster for large
%                   data sets)
%
%           NOTE: IF you want to pre-calculate the correction factor matrix
%           (which is recommended because it saves time), then run the
%           function:
%           [corrFacMatrix] = makeCorrFactorMatrix(imsiz, dist, samplesize,
%           mask);
%
%           IF nargin == 5              the provided corrFacMat is used
%           IF corrFacMat is an integer corrFacMat is calculated using
%                                       this integer value as samplesize
%           IF nargin == 4              Ripley's edge correction is used
%           IF nargin == 3 (or ==nan)   dist is set as default [1:rs]
%
% OUTPUT
%           kr:     Ripley's K-function
%           lr:     Besag's L-function = sqrt(K(r))-r
%           pcr:    pair-correlation function
%           glr:     Getis' L-function = sqrt(K(r))
%           for all, every column contains the kr/lr function for one
%           frame of the mpm-file, the row values correspond to the
%           specified distances
%
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
%
% last modified
% DATE:     29-Jan-2008 (last update)
%           1-July-2009
%
%

% create vector containing x- and y-image size
imsizex = imsiz(1);
imsizey = imsiz(2);

% if no distance vector is chosen explicitly or if dist is nan, the
% distance vector is set as default to 1:rs, where rs is the half diagonal
% of the image (this is the standard in the literature)
rs = round(sqrt(imsizex^2+imsizey^2));
if nargin<4
    distvec = 1:rs;
    nr = rs;
elseif isnan(dist)
    distvec = 1:rs;
    nr = rs;
else
    distvec = dist;
    nr = length(distvec);
end

% if corrFacMat is specified but consists only of an integer value,
% calculate the matrix using that value as sample, else set corrFacMat to
% empty
if (nargin<5)
    corrFacMat = [];
else
    if length(corrFacMat)==1
        sample = corrFacMat;
        corrFacMat = makeCorrFactorMatrix(imsiz,sample);
    end
end


if nargin<7
    sparse = false;
end

%determine size of mpm-file
[nx1,ny1]=size(mpm1);
%number of frames
numframes1 = round(ny1/2);
%determine size of mpm-file
[nx2,ny2]=size(mpm2);
%number of frames
numframes2 = round(ny2/2);

% function requires both mpms to have the same number of frames
if numframes1~=numframes2
    error('number of frames in 2 mpms doesn''t match');
else
    numf = numframes1;
end

%initialize results matrix pvr; x-dimension equals the employed number of
%values for the circle radius, y-dimension equals number of planes of the
%input mpm-file
kr  = zeros(nr,numf);
lr  = kr;
pcr = lr;
glr  = lr;

% loop over all frames
for i=1:numf
    % relevant nonzero x,y-positions in mpm1
    pos1    = find(mpm1(:,i*2)>0);
    cmpm1   = mpm1(pos1,2*i-1:2*i);
    np1     = length(pos1);
    
    % relevant nonzero x,y-positions in mpm2
    pos2    = find(mpm2(:,i*2)>0);
    cmpm2   = mpm2(pos2,2*i-1:2*i);
    np2     = length(pos2);
    
    %fprintf(' frame %04d',i);
    
    % if there are any relevant points in this frame - at least one point
    % each for cross-correlation (different matrices mpm1/mpm2), at least
    % two points for auto-corr (mpm1 == mpm2)
    if ( min(np1,np2)>0 ) && ( (np1+np2)>2 )
        
        % output pr: # of points as a function of distance
        [pr,nump] = pointsincircleCross(cmpm1,cmpm2,imsiz,distvec,corrFacMat,sparse);
        
        % normalized pr - normalize by total point density
        totaldensity = (nump)/(imsizex*imsizey);
        if nargin>5 && ~isempty(normArea)
            totaldensity = nump/normArea;
        end
        
        prnorm = (pr/pi)*(1/totaldensity);
        %prnorm = pr;
        kr(:,i)     = prnorm;
        lr(:,i)     = sqrt(prnorm) - distvec;
        pcr(:,i)    = convertLR2PCF(lr(:,i),distvec);
        glr(:,i) = sqrt(prnorm);
    else
        kr(:,i)     = nan*dist;
        lr(:,i)     = nan*dist;
        pcr(:,i)    = nan*dist;
        glr(:,i)    = nan*dist;
    end
    
    %fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for i-loop

%fprintf('\n');

end % of function





%=========================================================================
%=========================================================================
%=========================================================================
%====================       SUBFUNCTIONS    ==============================
%=========================================================================
%=========================================================================
%=========================================================================




function [npvr,nump]=pointsincircleCross(m1,m2,ms,dist,corrMat,sparse)
% pointsincircle calculates the average number of points in a circle around
% a given point as a function of the circle radius (averaged over all points
% and normalized by total point density); this function is called Ripley's
% K-function in statistics, and is an indication of the amount of clustering
% in the point distribution
%
% SYNOPSIS   [m2,num]=pointsincircle(m1,m2,ms,dist,corrMat);
%
% INPUT     m1:     matrix of size (n1 x 2) containing the (x,y)-coordinates
%                   of n1 points; these points are considered the CHILDREN
%           m2:     matrix of size (n2 x 2) containing the (x,y)-coordinates
%                   of n2 points; these points are considered the PARENTS
%           ms:     vector containing the parameters [imsizex imsizey] (the
%                   x-size and y-size of the image)
%           dist:   distance vector
%           corrMat:   correctionFactor matrix; if this matrix is empty or
%                   nargin<5, then the simple Ripley's edge correction is
%                   used to calculate the correction factor
%           sparse: true to measure distances using graph based algorithm
%                   that creates sparse matrix (might be faster for large
%                   data sets)
%
%
% OUTPUT    npvr:   vector containing the number of points in a circle
%                   around each point, for an increasing radius;
%                   radius default values are 1,2,3,....,min(ms)
%                   function is averaged over all objects in the image
%           nump:   number of points
%
% Dinah Loerke, Jan 29th, 2008

nump1 = size(m1,1);
msx=ms(1);
msy=ms(2);

% in the following, the matrix mdist will contain the distance of all
% points in m1 from all points in m2, where m1 are the children and m2 are
% the parents; because of the way the function below is set up in terms of
% rows-columns, the order in the DistanceMatrix function needs to be
% parent-child

%Distance Matrix is faster for smaller matrices, but
%createSparseDistanceMatrix can better deal with large matrices
if sparse
    mdist = createSparseDistanceMatrix([m2(:,2) m2(:,1)],[m1(:,2) m1(:,1)],max(dist(:)),0);
else
    [mdist]=DistanceMatrix(m2,m1);
end

% NOTE: In the old version of the Ripley, since it was designed for
% self-correlation, only distances > 0 were considered in the distance
% histogram below. In this version, since it allows cross-correlation,
% zero-distances have to be included for DIFFERENT mpms, but should be
% excluded for IDENTICAL mpms

% monitor progress
%fprintf(' progress %02d',0);

if isempty(m1) || isempty(m2) || isempty(mdist)
    npvr = nan*dist;
    nump = 0;
else
    
    %allocate space
    histmat = nan(length(mdist(:,1)),length(dist)+1);
    corrmat = nan(length(mdist(:,1)),length(dist));
    for i=1:length(mdist(:,1))
        
        % matrix histmat contains the distance histogram for each point, where
        % every row represents the distance histogram for the cell at that
        % position
        histmat(i,:) = histc(nonzeros(mdist(i,:)),[0 dist]);
        
        % determine the circumference correction vector for this point; in the
        % matrix corrmat, the columns represent the points, and the rows
        % represent the entries for the radius vector dist
        % IF corrMat exists, use these values directly, else calculate corrmat
        % in situ with Ripley's edge correction
        if ~isempty(corrMat)
            cpx = round(m2(i,1)); cpy = round(m2(i,2));
            corrmat(i,:)= corrMat(max(cpx,1),max(cpy,1),:);
        else
            % x- and y-distances of this center point (the parent point) from the
            % nearest edge
            ex = min(m2(i,1),1+msx-m2(i,1));
            ey = min(m2(i,2),1+msy-m2(i,2));
            corrmat(i,:)= circumferenceCorrFactor(ex,ey,dist,msx,msy);
        end
        
        %         % update iter
        %         iter = round( 100*(i/length(mdist(:,1))) );
        %         %if iter<100, fprintf('\b\b%02d',iter); end
        
    end
    
    
    
    % multiply the matrices histmat and the correction vector
    histmat(:,length(dist)+1)=[];
    pointsmat = histmat./corrmat;
    
    % for every point, the (number-of-points within radius r)-function is the
    % cumulative sum of the values in the corresponding row
    pointsincirclemat = cumsum(pointsmat,2);
    
    % average for this frame over all existing points
    npvr = nanmean(pointsincirclemat,1);
    nump = nump1;
    
end % of if there are any points



end % of function


function [m2]=DistanceMatrix(c1,c2)
%this subfunction makes a neighbour-distance matrix for input matrix c1
%(n1 x 2 points) and c2
%output: m2 (n1 x n1) matrix containing the distances of each point in c1
%from each point in c2

np1=size(c1,1);
np2=size(c2,1);

m2=zeros(np1,np2);

for k = 1:np1
    for n = 1:np2
        d = sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end % of for


end % of subfunction


%% =======================================================================
function pcfunc = convertLR2PCF(lr,dvec)
% convert L-function to pair correlation function

pcfunc = lr;
[dlx,dly] = size(lr);

if dlx==length(dvec)
    nf = dly;
else
    nf = dlx;
    lr = lr';
end

% area of circles corresponding to radii in dvec
area        = dvec.^2;
% area of radius increments (central circle and rings)
area_diff   = area;
area_diff(2:length(area_diff)) = diff(area);

% matrix of areas
amat    = repmat(area_diff',1,dly);
dmat    = repmat(dvec',1,dly);

for n=1:nf
    kr          = (lr(:,n)+dmat).^2;
    kr_diff     =  kr;
    kr_diff(2:length(dvec),:) = diff(kr,1);
    
    pcfunc(:,n) = kr_diff./amat;
end % of for

end % of subfunction
np2=size(c2,1);

m2=zeros(np1,np2);

for k = 1:np1
    for n = 1:np2
        d = sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end % of for


end % of subfunction


%% =======================================================================
function pcfunc = convertLR2PCF(lr,dvec)
% convert L-function to pair correlation function

pcfunc = lr;
[dlx,dly] = size(lr);

if dlx==length(dvec)
    nf = dly;
else
    nf = dlx;
    lr = lr';
end

% area of circles corresponding to radii in dvec
area        = dvec.^2;
% area of radius increments (central circle and rings)
area_diff   = area;
area_diff(2:length(area_diff)) = diff(area);

% matrix of areas
amat    = repmat(area_diff',1,dly);
dmat    = repmat(dvec',1,dly);

for n=1:nf
    kr          = (lr(:,n)+dmat).^2;
    kr_diff     =  kr;
    kr_diff(2:length(dvec),:) = diff(kr,1);
    
    pcfunc(:,n) = kr_diff./amat;
end % of for

end % of subfunction
