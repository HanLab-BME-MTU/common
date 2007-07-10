function [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
    nnDistFeatures,errFlag] = linkFeaturesKalman(movieInfo,costMatParam,...
    filterInfoPrev,kalmanInitParam,useLocalDensity,nnWindow,probDim)
%LINKFEATURESKALMAN links features between consecutive frames using LAP and directed motion models
%
%SYNOPSIS [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
%    nnDistFeatures,errFlag] = linkFeaturesKalman(movieInfo,costMatParam,...
%    filterInfoPrev,kalmanInitParam,useLocalDensity,nnWindow,probDim)
%
%INPUT  movieInfo         : Array of size equal to the number of frames
%                           in a movie, containing the fields:
%             .xCoord          : [Image coordinate system] x-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%             .yCoord          : [Image coordinate system] y-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%             .zCoord          : [Image coordinate system] z-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%                                Optional. Skipped if problem is 2D. Default: zeros.
%             .amp             : Amplitudes of PSFs fitting detected features. 
%                                1st column for values and 2nd column 
%                                for standard deviations.
%             .num             : Number if features. Optional. Calculated 
%                                if not supplied.  
%             .allCoord        : x,dx,y,dy,[z,dz] of features collected in one
%                                matrix. Optional. Calculated if not
%                                supplied.
%             .nnDist          : Distance from each feature to its nearest
%                                neighbor. Optional. Calculated if not
%                                supplied.
%       costMatParam      : Parameters needed for cost matrix calculation. 
%                           Structure with fields specified by particular
%                           cost matrices used.
%       filterInfoPrev    : Structure array with number of entries equal to
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%                           Optional. Enter [] or nothing if not to be used.
%       kalmanInitParam   : Structure with fields containing variables
%                           used in Kalman filter initialization. See
%                           particular initialization function for fields.
%                           Optional. Enter [] or nothing if not to be used.
%       useLocalDensity   : Logical variable indicating whether to use
%                           local density in search radius estimation.
%                           Optional. Default: 0.
%       nnWindow          : Time window to be used in estimating the
%                           nearest neighbor distance of a feature in a
%                           track. Needed even if useLocalDensity = 0. If
%                           not input, default is 1. 
%       probDim           : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%                           Optional. If not input, dimensionality will be
%                           derived from movieInfo.
%
%OUTPUT trackedFeatureIndx: Connectivity matrix of features between time points.
%                           Rows indicate continuous tracks, while columns 
%                           indicate frames. A track that ends before the
%                           last frame is followed by zeros, and a track
%                           that starts at a time after the first frame
%                           is preceded by zeros. 
%       trackedFeatureInfo: The positions and amplitudes of the tracked
%                           features. Number of rows = number of tracks, 
%                           while number of columns = 8*number of time 
%                           points. Each row consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%       kalmanFilterInfo  : Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       nnDistFeatures    : Matrix indicating the nearest neighbor
%                           distances of features linked together within
%                           tracks.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must not be empty.
%
%Khuloud Jaqaman, March 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackedFeatureIndx = [];
trackedFeatureInfo = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--linkFeaturesKalman: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether z-coordinates were input, making problem potentially 3D
if isfield(movieInfo,'zCoord')
    probDimT = 3;
else    
    probDimT = 2;
end

%check whether a priori Kalman filter information is given
if nargin < 3 || isempty(filterInfoPrev)
    filterInfoPrev = [];
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%check whether additional parameters for Kalman filter initialization are
%given
if nargin < 4 || isempty(kalmanInitParam)
    kalmanInitParam = [];
end

%check whether the variable useLocalDensity was input
if nargin < 5 || isempty(useLocalDensity)
    useLocalDensity = 0; %default is zero
end

%check whether the variable nnWindow was input
if nargin < 6 || isempty(nnWindow)
    nnWindow = 1;
end

%assign problem dimensionality if not input
if nargin < 7 || isempty(probDim)
    probDim = probDimT;
else
    if probDim == 3 && probDimT == 2
        disp('--linkFeaturesKalman: Inconsistency in input. Problem 3D but no z-coordinates.');
        errFlag = 1;
    end
end

%exit if there are problems with input
if errFlag
    disp('--linkFeaturesKalman: Please fix input parameters.');
    return
end

%get number of frames in movie
numFrames = length(movieInfo);

%calculate number of features in each frame if not supplied
if ~isfield(movieInfo,'num')
    for iFrame = 1 : numFrames
        movieInfo(iFrame).num = size(movieInfo(iFrame).xCoord,1);
    end
end

%collect coordinates in one matrix if allCoord not supplied
if ~isfield(movieInfo,'allCoord')
    if probDim == 2
        movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
            movieInfo(iFrame).yCoord];
    else
        movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
            movieInfo(iFrame).yCoord movieInfo(iFrame).zCoord];
    end
end

%calculate nearest neighbor distances if not supplied
if ~isfield(movieInfo,'nnDist')
    
    %calculate distance between each feature and its nearest neighbor
    %for all frames
    for iFrame = 1 : numFrames

        switch movieInfo(iFrame).num

            case 0 %if there are no features

                %there are no nearest neighbor distances
                nnDist = zeros(0,1);

            case 1 %if there is only 1 feature

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000;

            otherwise %if there are more than 1 feature

                %compute distance matrix
                nnDist = createDistanceMatrix(movieInfo(iFrame).allCoord(1:2:end),...
                    movieInfo(iFrame).allCoord(1:2:end));

                %sort distance matrix and find nearest neighbor distance
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,2);

        end

        %store nearest neighbor distance
        movieInfo(iFrame).nnDist = nnDist;

    end %(for iFrame = 1 : numFrames)

end %(if ~isfield(movieInfo,'nnDist'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reserve memory for kalmanFilterInfo - generalize to 3D - maybe make
%separate function
for iFrame = numFrames : -1 : 1
    numFeatures(iFrame) = movieInfo(iFrame).num;
end
kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim);

%fill the feature indices in 1st frame in the connectivity matrix
trackedFeatureIndx = (1:movieInfo(1).num)';

%fill the nearest neighbor distances of features in first frame, which are
%the first approximation of the nearest neighbor distances of tracks
nnDistFeatures = movieInfo(1).nnDist;
nnDistTracks = movieInfo(1).nnDist;

%initialize Kalman filter for features in 1st frame
if usePriorInfo %use a priori information if available
    kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
    kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
    kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
else
    [filterInit,errFlag] = kalmanInitLinearMotion(movieInfo(1),probDim,...
        kalmanInitParam);
    kalmanFilterInfo(1).stateVec = filterInit.stateVec;
    kalmanFilterInfo(1).stateCov = filterInit.stateCov;
    kalmanFilterInfo(1).noiseVar = filterInit.noiseVar;
end

%go over all frames
for iFrame = 1 : numFrames-1

    %get number of features
    numFeaturesFrame1 = movieInfo(iFrame).num; %in 1st frame
    numFeaturesFrame2 = movieInfo(iFrame+1).num; %in 2nd frame

    if numFeaturesFrame1 ~= 0 %if there are features in 1st frame

        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame

            %calculate cost matrix
            %function also outputs Kalman filter predictions for feature positions
            %and velocities in 2nd frame based on possible motion models
            [costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker] ...
                = costMatLinearMotionLink(movieInfo(iFrame:iFrame+1),...
                kalmanFilterInfo(iFrame),costMatParam,useLocalDensity,...
                nnDistTracks(1:numFeaturesFrame1),probDim);

            if any(costMat(:)~=nonlinkMarker) %if there are potential links

                %link features based on cost matrix, allowing for birth and death
                [link12,link21] = lap(costMat,nonlinkMarker,0,1);

                %get indices of features in 2nd frame that are connected to features in 1st frame
                indx2C = find(link21(1:numFeaturesFrame2)<=numFeaturesFrame1);

                %get indices of corresponding features in 1st frame
                indx1C = link21(indx2C);

                %find existing tracks that are not connected to features in 2nd frame
                indx1U = ones(size(trackedFeatureIndx,1),1);
                indx1U(indx1C) = 0;
                indx1U = find(indx1U);

                %assign space for new connectivity matrix,
                %for matrix of nearest neighbors of connected features,
                %and for vector storing track nearest neighbor distances
                tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2-length(indx2C),iFrame+1);
                tmpNN = NaN*ones(size(tmp));
                tmpNN2 = NaN*ones(size(tmp,1),1);

                %fill in the feature numbers in 2nd frame and their nearest
                %neighbor distances
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
                tmpNN(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

                %shuffle existing tracks to get the correct connectivity with 2nd frame
                %also shuffle nearest neighbor distances
                tmp(indx2C,1:iFrame) = trackedFeatureIndx(indx1C,:);
                tmpNN(indx2C,1:iFrame) = nnDistFeatures(indx1C,:);
                tmpNN2(indx2C) = nnDistTracks(indx1C);

                %add rows of tracks that are not connected to points in 2nd frame
                %also add nearest neighbor distances
                tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx(indx1U,:);
                tmpNN(numFeaturesFrame2+1:end,1:iFrame) = nnDistFeatures(indx1U,:);
                tmpNN2(numFeaturesFrame2+1:end) = nnDistTracks(indx1U);

                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;
                
                %update the nearest neighbor distances of connected
                %features and the track nearest neighbor vector
                nnDistFeatures = tmpNN;
                nnDistTracks = tmpNN2;
                tmpNN = max(1,iFrame-nnWindow);
                nnDistTracks(1:numFeaturesFrame2) = min(nnDistFeatures...
                    (1:numFeaturesFrame2,tmpNN:end),[],2);

                %use the Kalman gain from linking to get better estimates
                %of the state vector and its covariance matrix in 2nd frame
                %as well as state noise and its variance
                if usePriorInfo %if prior information is supplied
                    [kalmanFilterInfo,errFlag] = kalmanGainLinearMotion(...
                        trackedFeatureIndx(1:numFeaturesFrame2,:),...
                        movieInfo(iFrame+1),kalmanFilterInfoTmp,...
                        propagationScheme,kalmanFilterInfo,probDim,...
                        filterInfoPrev(iFrame+1),kalmanInitParam);
                else %if no prior information is supplied
                    [kalmanFilterInfo,errFlag] = kalmanGainLinearMotion(...
                        trackedFeatureIndx(1:numFeaturesFrame2,:),...
                        movieInfo(iFrame+1),kalmanFilterInfoTmp,...
                        propagationScheme,kalmanFilterInfo,probDim,[],...
                        kalmanInitParam);
                end
                
            else %if there are no potential links

                %assign space for new connectivity matrix,
                %for matrix of nearest neighbors of connected features,
                %and for vector storing track nearest neighbor distances
                tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2,iFrame+1);
                tmpNN = NaN*ones(size(tmp));
                tmpNN2 = NaN*ones(size(tmp,1),1);

                %fill in the feature numbers in 2nd frame and their nearest
                %neighbor distances
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
                tmpNN(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

                %fill in the tracks upto 1st frame and their nearest
                %neighbor distances
                tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx;
                tmpNN(numFeaturesFrame2+1:end,1:iFrame) = nnDistFeatures;
                tmpNN2(numFeaturesFrame2+1:end) = nnDistTracks;

                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;

                %update the nearest neighbor distances of connected
                %features and the track nearest neighbor vector
                nnDistFeatures = tmpNN;
                nnDistTracks = tmpNN2;
                nnDistTracks(1:numFeaturesFrame2) = movieInfo(iFrame+1).nnDist;

                %initialize Kalman filter for features in 2nd frame
                if usePriorInfo %use a priori information if available
                    kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                    kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                else
                    [filterInit,errFlag] = kalmanInitLinearMotion(...
                        movieInfo(iFrame+1),probDim,kalmanInitParam);
                    kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                    kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                end

            end %(if any(costMat(:)~=nonlinkMarker))

        else %if there are no features in 2nd frame

            %add a column of zeros for the 2nd frame in the track
            %connectivity matrix
            trackedFeatureIndx = [trackedFeatureIndx zeros(size(trackedFeatureIndx,1),1)];
            
            %add a column of NaNs for the 2nd frame in the nearest neighbor
            %matrix
            nnDistFeatures = [nnDistFeatures NaN*ones(size(nnDistFeatures,1),1)];

        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    else %if there are no feature in 1st frame

        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame

            %assign space for new connectivity matrix,
            %for matrix of nearest neighbors of connected features,
            %and for vector storing track nearest neighbor distances
            tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2,iFrame+1);
            tmpNN = NaN*ones(size(tmp));
            tmpNN2 = NaN*ones(size(tmp,1),1);

            %fill in the feature numbers in 2nd frame and their nearest
            %neighbor distances
            tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
            tmpNN(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

            %fill in the tracks upto 1st frame and their nearest
            %neighbor distances
            tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx;
            tmpNN(numFeaturesFrame2+1:end,1:iFrame) = nnDistFeatures;
            tmpNN2(numFeaturesFrame2+1:end) = nnDistTracks;

            %update the connectivity matrix "trackedFeatureIndx"
            trackedFeatureIndx = tmp;

            %update the nearest neighbor distances of connected
            %features and the track nearest neighbor vector
            nnDistFeatures = tmpNN;
            nnDistTracks = tmpNN2;
            nnDistTracks(1:numFeaturesFrame2) = movieInfo(iFrame+1).nnDist;

            %initialize Kalman filter for features in 2nd frame
            if usePriorInfo %use a priori information if available
                kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
            else
                [filterInit,errFlag] = kalmanInitLinearMotion(...
                    movieInfo(iFrame+1),probDim,kalmanInitParam);
                kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
            end

        else %if there are no features in 2nd frame

            %add a column of zeros for the 2nd frame in the track
            %connectivity matrix
            trackedFeatureIndx = [trackedFeatureIndx zeros(size(trackedFeatureIndx,1),1)];
            
            %add a column of NaNs for the 2nd frame in the nearest neighbor
            %matrix
            nnDistFeatures = [nnDistFeatures NaN*ones(size(nnDistFeatures,1),1)];

        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    end %(if numFeaturesFrame1 ~= 0 ... else ...)

end %(for iFrame=1:numFrames-1)

%get total number of tracks
numTracks = size(trackedFeatureIndx,1);

%find the frame where each track begins and then sort the vector
frameStart = zeros(numTracks,1);
for i=1:numTracks
    frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
end
[frameStart,indx] = sort(frameStart);

%rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the 
%same time in descending order from longest to shortest.
trackedFeatureIndx = trackedFeatureIndx(indx,:);

%also re-arrange the matrix indicating nearest neighbor distances
nnDistFeatures = nnDistFeatures(indx,:);

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system
%since this code is restricted to 2D data, z=dz=0

%reserve space for matrix
trackedFeatureInfo = NaN*ones(size(trackedFeatureIndx,1),8*numFrames);

%for now, the matrix always has space for z, hence I have to treat 2D and
%3D differently.
%at some point I can change the postprocessing functions to handle 2D and
%3D differently, and so I can save them differently, in which case I won't
%need the if-else here.
if probDim ==2

    %go over all frames
    for iFrame = 1 : numFrames

        %find rows that have a feature index
        indx1 = find(trackedFeatureIndx(:,iFrame)~=0);

        %if there are detected features in this frame
        if ~isempty(indx1)

            %find the feature index
            indx2 = trackedFeatureIndx(indx1,iFrame);

            %store its information
            trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
                [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
                zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,1) ...
                movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
                zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,2)];

        end

    end

else

    %go over all frames
    for iFrame = 1 : numFrames

        %find rows that have a feature index
        indx1 = find(trackedFeatureIndx(:,iFrame)~=0);

        %if there are detected features in this frame
        if ~isempty(indx1)

            %find the feature index
            indx2 = trackedFeatureIndx(indx1,iFrame);

            %store its information
            trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
                [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
                movieInfo(iFrame).amp(indx2,1) ...
                movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
                movieInfo(iFrame).amp(indx2,2)];

        end

    end

end

%%%%% ~~ the end ~~ %%%%%
