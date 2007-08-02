function [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal)

%get number of tracks
numTracks = length(tracksFinal);

%get number of time points
tmp = vertcat(tracksFinal.seqOfEvents);
numTimePoints = max(tmp(:,1));

%reserve memory for matrix of tracks
trackedFeatureInfo = NaN(numTracks,8*numTimePoints);

%reserve memory for matrix of feature indices
trackedFeatureIndx = zeros(numTracks,numTimePoints);

%put tracks in matrix
for iTrack = 1 : numTracks
    startTime = tracksFinal(iTrack).seqOfEvents(1,1);
    endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
    trackedFeatureInfo(iTrack,8*(startTime-1)+1:8*endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG;
    trackedFeatureIndx(iTrack,startTime:endTime) = ...
        tracksFinal(iTrack).tracksFeatIndxCG;
end
