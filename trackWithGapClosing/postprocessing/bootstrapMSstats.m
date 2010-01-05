function msBootstrapRes = bootstrapMSstats(tracksAll,minTrackLen,probDim,...
    diffAnalysisAll,numRep)

%remove tracks shorter than minTrackLen
criteria.lifeTime.min = minTrackLen;
indxKeep = chooseTracks(tracksAll,criteria);
tracksAll = tracksAll(indxKeep);
diffAnalysisAll = diffAnalysisAll(indxKeep);

%get number of available tracks
numTracks = length(tracksAll);

%reserve memory for results
msBootstrapRes = NaN(numRep,27);

%repeat numRep times
for iRep = 1 : numRep

    %select a random subset, with replacement
    if iRep ~= 1
        randIndx = randsample((1:numTracks),numTracks,'true');
    else
        randIndx = 1:numTracks;
    end
    tracksSample = tracksAll(randIndx);
    diffAnalysisSample = diffAnalysisAll(randIndx);

    %analyze MS stats
    [statsGeneral,statsPerCat] = calcStatsMS(tracksSample,minTrackLen,...
        probDim,diffAnalysisSample);

    %extract results
    probMS = statsGeneral(4:5); %overall merging and splitting probabilities
    probType = statsPerCat(:,2)'; %overall type probabilities
    tmp = statsPerCat(:,3:4)';
    prob1 = tmp(:)'; %conditional probabilities as calculated in approach 1
    tmp = statsPerCat(:,5:6)';
    prob2 = tmp(:)'; %conditional probabilities as calculated in approach 2
% % %     tmp1 = statsPerCat(:,7:10);
% % %     tmp1 = tmp1(:)';
% % %     tmp2 = statsPerCat(:,11:14);
% % %     tmp2 = tmp2(:)';
% % %     tmp1 = [tmp1; tmp2];
% % %     tmp1 = tmp1(:,[1 5 9 13 6 10 14 11 15 16]);
% % %     prob3 = tmp1(:)'; %conditional probabilities as calculated in approach 3

    %store results
    msBootstrapRes(iRep,:) = [probMS probType prob1 prob2]; % prob3];

end

