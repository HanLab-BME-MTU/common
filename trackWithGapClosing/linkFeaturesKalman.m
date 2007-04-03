function [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,errFlag] ...
    = linkFeaturesKalman(movieInfo,costMatParam,filterInfoPrev)
%LINKFEATURESKALMAN links features between consecutive frames using LAP and directed motion models
%
%SYNOPSIS [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,errFlag] ...
%    = linkFeaturesKalman(movieInfo,linkParam,filterInfoPrev)
%
%INPUT  movieInfo         : Array of size equal to the number of time points
%                           in a movie, containing the fields:
%             .xCoord          : Image coordinate system x-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%             .yCoord          : Image coordinate system y-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%             .amp             : Amplitudes of PSFs fitting detected features. 
%                                1st column for values and 2nd column 
%                                for standard deviations.
%       costMatParam      : Parameters needed for cost matrix calculation. 
%                           Structure with fields specified by particular
%                           cost matrix.
%       filterInfoPrev    : Structure array with number of entries equal to
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%                           Optional input. Enter [] or nothing if not to be used.
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
%
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS Code restricted to 2D motion. 
%        Algorithm can handle cases where some frames do not have any
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

%check whether a priori Kalman filter information is given
if nargin < 3 || isempty(filterInfoPrev)
    filterInfoPrev = [];
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of frames in movie
numFrames = length(movieInfo);

%get number of features in each frame
for iFrame = 1 : numFrames
    movieInfo(iFrame).num = size(movieInfo(iFrame).xCoord,1);
end

%reserve memory for kalmanFilterInfo
for iFrame = 1 : numFrames
    numFeatures = movieInfo(iFrame).num;
    kalmanFilterInfo(iFrame) = struct('stateVec',zeros(numFeatures,4),...
        'stateCov',zeros(4,4,numFeatures),'noiseVar',zeros(4,4,numFeatures),...
        'stateNoise',zeros(numFeatures,4),'scheme',zeros(numFeatures,2));
end

%fill the feature numbers in 1st frame in the connectivity matrix
trackedFeatureIndx = (1:movieInfo(1).num)';

%initialize Kalman filter for features in 1st frame
if usePriorInfo %use a priori information if available
    kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
    kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
    kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
else
    [filterInit,errFlag] = kalmanInitLinearMotion(movieInfo(1));
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
                = costMatLinearMotion(movieInfo(iFrame:iFrame+1),...
                kalmanFilterInfo(iFrame),costMatParam);

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

                %assign space for new matrix
                tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2-length(indx2C),iFrame+1);

                %fill in the feature numbers in 2nd frame
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';

                %shuffle existing tracks to get the correct connectivity with 2nd frame
                tmp(indx2C,1:iFrame) = trackedFeatureIndx(indx1C,:);

                %add rows of tracks that are not connected to points in 2nd frame
                tmp(max(numFeaturesFrame2,length(indx1C))+1:end,1:iFrame) ...
                    = trackedFeatureIndx(indx1U,:);

                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;

                %use the Kalman gain from linking to get better estimates
                %of the state vector and its covariance matrix in 2nd frame
                %as well as state noise and its variance
                if usePriorInfo %if prior information is supplied
                    [kalmanFilterInfo,errFlag] = kalmanGainLinearMotion(...
                        trackedFeatureIndx(1:numFeaturesFrame2,:),...
                        movieInfo(iFrame+1),kalmanFilterInfoTmp,...
                        propagationScheme,kalmanFilterInfo,filterInfoPrev(iFrame+1));
                else %if no prior information is supplied
                    [kalmanFilterInfo,errFlag] = kalmanGainLinearMotion(...
                        trackedFeatureIndx(1:numFeaturesFrame2,:),...
                        movieInfo(iFrame+1),kalmanFilterInfoTmp,...
                        propagationScheme,kalmanFilterInfo);
                end
                
            else %if there are no potential links

                %assign space for new matrix
                tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2,iFrame+1);

                %fill in the feature numbers in 2nd frame
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';

                %fill in the tracks upto 1st frame
                tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx;

                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;

                %initialize Kalman filter for features in 2nd frame
                if usePriorInfo %use a priori information if available
                    kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                    kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                else
                    [filterInit,errFlag] = kalmanInitLinearMotion(movieInfo(iFrame+1));
                    kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                    kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                end

            end %(if any(costMat(:)~=nonlinkMarker))

        else %if there are no features in 2nd frame

            %add a column of zeros for the 2nd frame
            trackedFeatureIndx = [trackedFeatureIndx zeros(size(trackedFeatureIndx,1),1)];

        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    else %if there are no feature in 1st frame

        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame

            %assign space for new matrix
            tmp = zeros(size(trackedFeatureIndx,1)+numFeaturesFrame2,iFrame+1);

            %fill in the feature numbers in 2nd frame
            tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';

            %fill in the tracks upto 1st frame
            tmp(numFeaturesFrame2+1:end,1:iFrame) = trackedFeatureIndx;
            
            %update the connectivity matrix "trackedFeatureIndx"
            trackedFeatureIndx = tmp;

            %initialize Kalman filter for features in 2nd frame
                if usePriorInfo %use a priori information if available
                    kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                    kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                else
                    [filterInit,errFlag] = kalmanInitLinearMotion(movieInfo(iFrame+1));
                    kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                    kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                end

        else %if there are no features in 2nd frame

            %add a column of zeros for 2nd frame
            trackedFeatureIndx = [trackedFeatureIndx zeros(size(trackedFeatureIndx,1),1)];

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

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system
%since this code is restricted to 2D data, z=dz=0

%reserve space for matrix
trackedFeatureInfo = NaN*ones(size(trackedFeatureIndx,1),8*numFrames);

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
            [movieInfo(iFrame).xCoord(indx2,1) movieInfo(iFrame).yCoord(indx2,1)  ...
            zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,1) ...
            movieInfo(iFrame).xCoord(indx2,2) movieInfo(iFrame).yCoord(indx2,2) ...
            zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,2)];

    end

end


%%%%% ~~ the end ~~ %%%%%
