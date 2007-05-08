function [trackType,xVel,yVel,noiseStd,errFlag] = estimTrackTypeParamLM(...
    trackedFeatIndx,trackedFeatInfo,kalmanFilterInfo,trackConnect,lenForClassify)
%ESTIMTRACKTYPEPARAMLM ...
%
%SYNOPSIS [trackType,xVel,yVel,noiseStd,errFlag] = estimTrackTypeParamLM(...
%    trackedFeatIndx,trackedFeatInfo,kalmanFilterInfo,trackConnect,lenForClassify);
%
%INPUT  trackedFeatIndx : Connectivity matrix of features between time
%                         points from initial linking. Rows indicate tracks, while columns
%                         indicate frames. A track that ends before the
%                         last frame is followed by zeros, and a track
%                         that starts at a time after the first frame
%                         is preceded by zeros. 
%       trackedFeatInfo : The positions and amplitudes of the tracked
%                         features from initial linking. Number of rows = number of tracks,
%                         while number of columns = 8*number of time points.
%                         Each row consists of
%                         [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                         in image coordinate system (coordinates in
%                         pixels). NaN is used to indicate time points
%                         where the track does not exist.
%       kalmanFilterInfo: Kalman filter information as calculated in
%                         linkFeaturesKalman for the linear motion model.
%       trackConnect    : Matrix indicating connectivity between tracks from
%                         initial linking, for example due to gap closing.
%       lenForClassify  : Minimum length of a track to classify it as
%                         directed or Brownian.
%
%OUTPUT trackedFeatKalmanInfo: The velocities, random element variance and
%                              propagation scheme along the tracks. Rows
%                              indicate tracks. # columns = 8 * # time
%                              points. Each row consists of
%                              [vx1 dvx1 vy1 dvy1 epsX1 epsY1 varEps1 scheme1 vx2 dvx2 vy2 dvy2 epsX2 epsY2 varEps2 scheme2 ...]
%
%       errFlag              : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackType = [];
xVel = [];
yVel = [];
noiseStd = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('estimTrackTypeParamLM')
    disp('--estimTrackTypeParamLM: Incorrect number of input arguments!');
    return
end

%get number of tracks from initial linking and number of frames
[numTracksLink,numFrames] = size(trackedFeatIndx);

%get number of compound tracks in trackConnect
numTracksCG = size(trackConnect,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation of type and motion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reserve memory for output variables
trackType = NaN*ones(numTracksLink,1);
xVel = zeros(numTracksLink,1);
yVel = zeros(numTracksLink,1);
noiseStd = zeros(numTracksLink,1);

%get the start times, end times and lifetimes of all tracks
trackSEL = getTrackSEL(trackedFeatInfo);
trackStartTime = trackSEL(:,1);
trackEndTime   = trackSEL(:,2);
trackLifeTime  = trackSEL(:,3);

%assign the asymmetry parameter thresholds that indicate directed motion
%for different track lengths 
%90th percentile:
asymThresh = [[NaN NaN 5 2.8 2.2 1.9 1.7 1.6 1.5 1.5 1.45 1.4 1.4 1.4 1.4 1.4 ...
    1.4 1.35 1.35 1.3]'; 1.3*ones(numFrames-20,1)];
% % % %99th percentile:
% % % asymThresh = [[NaN NaN 10 5 3.7 3 2.8 2.7 2.6 2.5 2.4 2.3 2.2 2.2 2.1 2.1 ...
% % %     2.1 2.1 2.1 2.1]'; 2*ones(numFrames-20,1)];

%go over all compound tracks in trackConnect
for iTrack = 1 : numTracksCG

    %get the indices of track segments making this compound track
    segmentIndx = trackConnect(iTrack,:);
    segmentIndx = segmentIndx(segmentIndx ~= 0);
    numSegments = length(segmentIndx);
    
    %reserve memory for track asymmetry
    asymmetry = NaN*ones(numSegments,1); %for individual segments

    %go over all segments
    for i = 1 : numSegments

        %get index of current segment
        iSegment = segmentIndx(i);

        %if segments's lifetime is at least "lenForClassify" frames
        %(shorter tracks are not reliable) ...
        if trackLifeTime(iSegment) >= lenForClassify

            %get current track's x and y coordinates
            currentTrack = [trackedFeatInfo(iSegment,1:8:end)' ...
                trackedFeatInfo(iSegment,2:8:end)'];
            currentTrack = currentTrack(trackStartTime(iSegment):trackEndTime(iSegment),:);

            %evaluate the asymmetry parameter as defined in the Huet et al. BJ 2006 paper
            asymmetry(i) = asymDetermination(currentTrack);

        end %(if trackLifeTime(iSegment) >= lenForClassify)

    end %(for i = 1 : numSegments)

    %get maximum asymmetry among all segments
    [maxAsymmetry,indxMaxAsym] = max(asymmetry);
    indxMaxAsym = segmentIndx(indxMaxAsym);

    %if the maximum asymmetry among all segments is larger than threshold ...
    if maxAsymmetry > asymThresh(trackLifeTime(indxMaxAsym))

        %get coordinates of compound track
        xCoord = reshape(trackedFeatInfo(segmentIndx,1:8:end)',numSegments*numFrames,1);
        yCoord = reshape(trackedFeatInfo(segmentIndx,2:8:end)',numSegments*numFrames,1);
        currentTrack = [xCoord yCoord];
        currentTrack = currentTrack(~isnan(xCoord),:);

        %get number of observation in compound track
        overallTrackLength = size(currentTrack,1);

        %evaluate asymmetry of whole track
        overallAsymmetry = asymDetermination(currentTrack);

        %if asymmetry of overall track is larger than threshold, overall track motion is
        %directed (1). If smaller, overall track motion is Brownian (0).
        overallType = overallAsymmetry > asymThresh(overallTrackLength);
        
    else %if maxAsymmetry is not larger than threshold ...
        
        %assign type to zero
        overallType = 0;
        
    end

    %assign types and motion parameters to segments based on overallType
    switch overallType

        case 1 %if compound track is directed
            
            %assign the type of all of its segments to 1 (directed)
            trackType(segmentIndx) = 1;
            
            %assign to all segments the velocity and std of the track with
            %maximum asymmetry
            xVel(segmentIndx) = kalmanFilterInfo(trackEndTime(...
                indxMaxAsym)).stateVec(trackedFeatIndx(indxMaxAsym,...
                trackEndTime(indxMaxAsym)),3);
            yVel(segmentIndx) = kalmanFilterInfo(trackEndTime(...
                indxMaxAsym)).stateVec(trackedFeatIndx(indxMaxAsym,...
                trackEndTime(indxMaxAsym)),4);
            noiseStd(segmentIndx) = sqrt(kalmanFilterInfo(trackEndTime(...
                indxMaxAsym)).noiseVar(1,1,trackedFeatIndx(indxMaxAsym,...
                trackEndTime(indxMaxAsym))));

        otherwise
            
            switch isnan(maxAsymmetry)

                case 0 %if compound track is Brownian

                    %assign the type of all of its segments to 0 (Brownian)
                    trackType(segmentIndx) = 0;

                    %give all segments a velocity of zero
                    xVel(segmentIndx) = 0;
                    yVel(segmentIndx) = 0;

                    %find the segments which were originally Brownian
                    indxBrown = segmentIndx(~isnan(asymmetry));

                    %get their noise std
                    segNoiseStd = zeros(length(indxBrown),1);
                    for i = 1 : length(indxBrown)
                        segNoiseStd(i) = sqrt(kalmanFilterInfo(trackEndTime(...
                            indxBrown(i))).noiseVar(1,1,trackedFeatIndx(indxBrown(i),...
                            trackEndTime(indxBrown(i)))));
                    end

                    %assign the maximum noise std to all segments
                    noiseStd(segmentIndx) = max(segNoiseStd);

                otherwise %if all segments are undetermined

                    %assign the type of all of its segments to NaN (undetermined)
                    trackType(segmentIndx) = NaN;

                    %give all segments a velocity of zero
                    xVel(segmentIndx) = 0;
                    yVel(segmentIndx) = 0;
                    
                    %give all segments a noise std of 1 (this value will be
                    %overwritten in the actual calculation of average displacement
                    %and search radius
                    noiseStd(segmentIndx) = 1;

            end %(switch isnan(overallAsymmetry))

    end %(switch overallType)

end %(for iTrack = 1 : numTracksCG)
