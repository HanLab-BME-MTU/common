function plotTrackStartEnd(tracksFinal)

tracksMat = convStruct2MatNoMS(tracksFinal);

xCoord = tracksMat(:,1:8:end);
yCoord = tracksMat(:,2:8:end);

trackSEL = getTrackSEL(tracksFinal);

numTracks = length(tracksFinal);

[xCoordStart,yCoordStart,xCoordEnd,yCoordEnd] = deal(NaN(numTracks,1));

for iTrack = 1 : numTracks
    xCoordStart(iTrack) = xCoord(iTrack,trackSEL(iTrack,1));
    yCoordStart(iTrack) = yCoord(iTrack,trackSEL(iTrack,1));
    xCoordEnd(iTrack) = xCoord(iTrack,trackSEL(iTrack,2));
    yCoordEnd(iTrack) = yCoord(iTrack,trackSEL(iTrack,2));
end

figure, hold on
for iTrack = 1 : numTracks
    plot([xCoordStart(iTrack) xCoordEnd(iTrack)],[yCoordStart(iTrack) yCoordEnd(iTrack)])
end
