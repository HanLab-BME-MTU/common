function [ combinedTrack ] = combine( obj )
%combine Combine multiple tracks together by uniquely indexing their
%segments

combinedTrack = TracksHandle;
combinedTrack.startFrame = min([obj.startFrame]);
combinedTrack.endFrame = max([obj.endFrame]);

segCount = cumsum([obj.numSegments]);
numSegments = segCount(end);
segCount = [0 segCount(1:end)] + 1;

tracksFeatIndxCG = zeros(numSegments,combinedTrack.lifetime);
tracksCoordAmpCG3D = NaN(numSegments,8,combinedTrack.lifetime);
seqOfEvents = obj.getSeqOfEventsMatrix;

for ii=1:numel(obj)
    tracksFeatIndxCG( ...
        segCount(ii):segCount(ii+1)-1, ...
        obj(ii).f-combinedTrack.startFrame+1) ...
        = obj(ii).tracksFeatIndxCG;
    
    tracksCoordAmpCG3D( ...
        segCount(ii):segCount(ii+1)-1, ...
        1:8, ...
        obj(ii).f-combinedTrack.startFrame+1) ...
        = obj(ii).tracksCoordAmpCG3D;
    
    combinedTrack.t(obj(ii).f) = obj(ii).t;
end

combinedTrack.tracksFeatIndxCG = tracksFeatIndxCG;
combinedTrack.tracksCoordAmpCG3D = tracksCoordAmpCG3D;
combinedTrack.seqOfEvents = seqOfEvents(:,1:4);


end

