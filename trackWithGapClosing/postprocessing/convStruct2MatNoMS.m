function trackedFeatureInfo = convStruct2MatNoMS(tracksFinal)

%get number of tracks
numTracks = length(tracksFinal);

%get number of time points
tmp = vertcat(tracksFinal.seqOfEvents);
numTimePoints = max(tmp(:,1));

%reserve memory for matrix of tracks
trackedFeatureInfo = NaN*ones(numTracks,8*numTimePoints);

%put tracks in matrix
for iTrack = 1 : numTracks
    startTime = tracksFinal(iTrack).seqOfEvents(1,1);
    endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
    trackedFeatureInfo(iTrack,8*(startTime-1)+1:8*endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG;
end
