%initialize temporary structures
tmpD = repmat(struct('field',[]),120,1);
tmpT = tmpD;

j = 0;
for startFrame = 1 : 400 : 48000
    endFrame = startFrame + 399;
    
    %get tracks for this time interval
    tracksFileName = ['tracks2Detection1_Frames' sprintf('%05i',startFrame) 'to' sprintf('%05i',endFrame) '.mat'];
    load(tracksFileName);
    
    %do diffusion analysis
    diffAnalysisRes = trackDiffusionAnalysis1(tracksFinal,1,2,1,[0.05 0.1],0,0);
    
    %save diffusion analysis of this time interval
    diffFileName = ['diffAnalysis1_Frames' sprintf('%05i',startFrame) 'to' sprintf('%05i',endFrame) '.mat'];
    save(diffFileName,'diffAnalysisRes');
    
    %shift the time stored in tracksFinal to reflect the real time
    numTracks = length(tracksFinal);
    for iTrack = 1 : numTracks
        tracksFinal(iTrack).seqOfEvents(:,1) = tracksFinal(iTrack).seqOfEvents(:,1) + startFrame - 1; %#ok<SAGROW>
    end
    
    %store tracks and diffusion analysis in temporary structures
    j = j + 1;
    tmpD(j).field = diffAnalysisRes;
    tmpT(j).field = tracksFinal;
    
end

%save combined diffusion analysis
diffAnalysisRes = vertcat(tmpD.field);
save('diffAnalysis1AllFrames','diffAnalysisRes');

%save combined tracks
tracksFinal = vertcat(tmpT.field);
save('tracks2Detection1AllFrames','tracksFinal');
