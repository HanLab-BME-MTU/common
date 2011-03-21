%initialize temporary structures
tmpD = repmat(struct('field',[]),17,1);
tmpT = tmpD;

j = 0;
for startFrame = 1 : 400 : 7200
    endFrame = startFrame + 399;
    
    %get tracks for this time interval
    tracksFileName = ['tracks8Detection2_Frames' sprintf('%04i',startFrame) 'to' sprintf('%04i',endFrame) '.mat'];
    load(tracksFileName);
    
    %do diffusion analysis
    diffAnalysisRes = trackDiffusionAnalysis1(tracksFinal,1,2,0,0.05,0,0);
    
    %save diffusion analysis of this time interval
    diffFileName = ['diffAnalysis_Frames' sprintf('%04i',startFrame) 'to' sprintf('%04i',endFrame) '.mat'];
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
save('diffAnalysisAllFrames','diffAnalysisRes');

%save combined tracks
tracksFinal = vertcat(tmpT.field);
save('tracks8Detection2AllFrames','tracksFinal');
