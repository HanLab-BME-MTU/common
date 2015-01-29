function tracks = normalizeTracks(tracks,movieInfo)
    if(~isfield(tracks,'seqOfEvents'))
        [tracks.seqOfEvents] = deal([1 1 1 NaN; length(movieInfo) 2 1 NaN]);
    end
    if(~isfield(tracks,'tracksCoordAmpCG'))
        tracksCoordAmpCG = getFeatFromIdx(tracks,movieInfo);
        [tracks.tracksCoordAmpCG] = deal(tracksCoordAmpCG{:});
    end
end
