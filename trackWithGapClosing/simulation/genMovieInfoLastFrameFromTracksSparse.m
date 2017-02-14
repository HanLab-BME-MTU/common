function compTracksAggregState = genMovieInfoLastFrameFromTracksSparse(tracksSim)
%GENMOVIEINFOFROMTRACKS generates a list of detected features per frame from supplied tracks
%
%SYNOPSIS movieInfoLastFrame = genMovieInfoFromTracks(tracksSim,percentMissing)
%
%INPUT  tracksSim     : simulate tracks (tracksSim).
%      
%OUTPUT movieInfoLastFrame: List of detected features in the last frame, in the format
%                  required for the input of trackWithGapClosing and
%                  trackCloseGapsKalman.
%
%                  .xCoord      : x-coordinates of detected features. 
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                  .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             
%                  .amp         : "Intensities" of detected features.
%                            1st column: values (ones if not available),
%                            2nd column: standard deviation (zeros if not
%                            available).
%
%
%Khuloud Jaqaman, October 2007
% Modified by Luciana de Oliveira, January 2017.

%call imput variables from structure
% tracksSim=tracksSim.tracksSim;

%get number of frames in movie
seqOfEvents = vertcat(tracksSim.seqOfEvents);
lastFrame = max(seqOfEvents(:,1));


%initialize variables
    xCoord= [];
    yCoord= [];
    amp   = [];

%transform each track in compTrack in the full matrix version
 
for     i=1:length(tracksSim); 
    tracksSimFull.tracksFeatIndxCG = full(tracksSim(i).tracksFeatIndxCG);
    tracksSimFull.tracksCoordAmpCG = full(tracksSim(i).tracksCoordAmpCG); 
    tracksSimFull.seqOfEvents = tracksSim(i).seqOfEvents; 
    [trackedFeatureInfo,~,~,~,~] = convStruct2MatIgnoreMS(tracksSimFull);
    
% take coordinates and intensity     
     xCoord = [xCoord; ...
                trackedFeatureInfo(:,(lastFrame-1)*8+1),trackedFeatureInfo(:,(lastFrame-1)*8+5)];
     yCoord = [yCoord; ...
                trackedFeatureInfo(:,(lastFrame-1)*8+2),trackedFeatureInfo(:,(lastFrame-1)*8+6)];
     amp    = [amp; ...
                trackedFeatureInfo(:,(lastFrame-1)*8+4),trackedFeatureInfo(:,(lastFrame-1)*8+8)];

end

% remove the tracks that are not present in the last frame (all row will be zero)
indx = find(any([xCoord yCoord amp],2));  
xCoord=xCoord(indx,:);
yCoord=yCoord(indx,:);
amp=amp(indx,:);
    
%store feature information in movieInfoLastFrame
    compTracksAggregState.xCoord = xCoord;
    compTracksAggregState.yCoord = yCoord;
    compTracksAggregState.amp = amp;
    

%%% ~~ the end ~~ %%%
