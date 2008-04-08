function [statsGeneral,statsPerCat,numEventsPerCat,msTimeInfo] = calcStatsMS(tracks,...
    minTrackLen,probDim,diffAnalysisRes)
%CALCSTATSMS calculates merge/split statistics for linear, Brownian and confined Brownian tracks
%
%SYNOPSIS [statsGeneral,statsPerCat,numEventsPerCat,msTimeInfo] = calcStatsMS(tracks,...
%    minTrackLen,probDim,diffAnalysisRes)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffAnalysis1). Optional. If not input, it will
%                    be calculated.
%
%OUTPUT statsGeneral: Row vector with entries: 
%                     (1) Number of features per frame.
%                     (2) Number of tracks with length >= minTrackLen.
%                     (3) Number of tracks segments belonging to the tracks
%                         with length >= minTrackLen.
%                     (4) Probability of a feature merging.
%                     (5) Probability of a feature splitting.
%       statsPerCat : 5-by-6 array. Rows correspond to: (1) linear, 
%                     (2) Brownian, (3) confined Brownian,
%                     (4) undetermined >= 5 frames, and 
%                     (5) undetermined < 5 frames.
%                     Columns corrrespond to:
%                     (1) Fraction of track segments in each category.
%                     (2) Probability of a feature belonging to a category.
%                     (3) Conditional probability of a feature
%                         merging IF in a category, where the category of
%                         the more dynamic feature is dominant.
%                     (4) Conditional probability of a feature
%                         splitting IF in a category, where the category of
%                         the more dynamic feature is dominant.
%                     (5) Conditional probability of a feature
%                         merging IF in a category, regardless of the other
%                         feature's category.
%                     (6) Conditional probability of a feature
%                         splitting IF in a category, regardless of the other
%                         feature's category.
% % %                     NO LONGER:
% % %                     (7-10) Conditional probability of a feature merging
% % %                            IF in a category AND the other feature is in
% % %                            the linear, Brownian, confined or
% % %                            undetermined category.
% % %                     (11-14) Conditional probability of a feature splitting
% % %                            IF in a category AND the other feature is in
% % %                            the linear, Brownian, confined or
% % %                            undetermined category.
%       numEventsPerCat: 5-by-4 array. Rows correspond to: (1) linear, 
%                     (2) Brownian, (3) confined Brownian,
%                     (4) undetermined >= 5 frames, and
%                     (5) undetermined < 5 frames.
%                     Columns corrrespond to:
%                     (1) Number of merges, where the category of the more 
%                         dynamic feature is dominant.
%                     (2) Number of splits, where the category of the more 
%                         dynamic feature is dominant.
%                     (3) Number of merges, regardless of the other
%                         feature's category.
%                     (4) Number of splits, regardless of the other
%                         feature's category.
% % %                     NO LONGER:
% % %                     (5-8) Number of merges when the other feature in
% % %                         linear, Brownian, confined or undetermined.
% % %                     (9-12) Number of splits when the other feature in
% % %                         linear, Brownian, confined or undetermined.
%       msTimeInfo : Structure with field 'brown' and 'linear' for Brownian
%                    and linear tracks. Each field is a structure
%                    containing the fields:
%           .numTracks      : 4 entries: # of tracks with merges and
%                             splits, # of tracks with only merges, # of
%                             tracks with only splits, and # of tracks
%                             without any merges or splits.
%           .timeMerge2Split: 2-column vector where first row is time
%                             between a merge and a consecutive split and
%                             second column is the weight of that
%                             observation.
%           .timeSplit2Merge: 2-column vector where first row is time
%                             between a split and a consecutive merge and
%                             second column is the weight of that
%                             observation.
%
%Khuloud Jaqaman, December 2007

%% input

if nargin < 1 || isempty(tracks)
    disp('calcStatsMS: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(diffAnalysisRes)
    [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,1,probDim,...
        1,[0.05 0.2],0);
    if errFlag
        return
    end
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);
diffAnalysisRes = diffAnalysisRes(indx);

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%put tracks in matrix format
[tracksMat,tracksIndxMat,trackStartRow,numSegments] = convStruct2MatIgnoreMS(tracks);

%get number of track segments
numTrackSegments = size(tracksMat,1);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% track segment types

%get track segment classification from diffusion analysis
trackSegmentClass = vertcat(diffAnalysisRes.classification);

%get track segment length
trackSegmentLength = getTrackSEL(tracksMat);
trackSegmentLength = trackSegmentLength(:,3);

%get indices of linear, Brownian, confined Brownian and undetermined tracks
indxLin    = find( trackSegmentClass(:,1) == 1 );
indxBrown  = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 2 );
indxConf   = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 1 );
indxUndet1 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
    & trackSegmentLength >= 5);
indxUndet2 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
    & trackSegmentLength < 5);

%store track segment type in an array
%-1 = undetermined & track length < 5, 0 = undetermined & track length >= 5, 
%1 = confined Brownian, 2 = Brownian, 3 = linear
trackSegmentType = zeros(numTrackSegments,1);
trackSegmentType(indxUndet2) = -1;
trackSegmentType(indxConf)   = 1;
trackSegmentType(indxBrown)  = 2;
trackSegmentType(indxLin)    = 3;

%calculate number of track segments per category
numSegmentsLin    = length(indxLin);
numSegmentsBrown  = length(indxBrown);
numSegmentsConf   = length(indxConf);
numSegmentsUndet1 = length(indxUndet1);
numSegmentsUndet2 = length(indxUndet2);

%calculate fraction of track segments falling in each category
fracSegmentsLin    = numSegmentsLin    / numTrackSegments;
fracSegmentsBrown  = numSegmentsBrown  / numTrackSegments;
fracSegmentsConf   = numSegmentsConf   / numTrackSegments;
fracSegmentsUndet1 = numSegmentsUndet1 / numTrackSegments;
fracSegmentsUndet2 = numSegmentsUndet2 / numTrackSegments;

%calculate number of features in each category
numFeatLin    = length(find(tracksIndxMat(indxLin,:)));
numFeatBrown  = length(find(tracksIndxMat(indxBrown,:)));
numFeatConf   = length(find(tracksIndxMat(indxConf,:)));
numFeatUndet1 = length(find(tracksIndxMat(indxUndet1,:)));
numFeatUndet2 = length(find(tracksIndxMat(indxUndet2,:)));

%get fraction of features in each category - this is the probability of a
%feature undergoing a certain motion type
probFeatLin    = numFeatLin    / numFeatTot;
probFeatBrown  = numFeatBrown  / numFeatTot;
probFeatConf   = numFeatConf   / numFeatTot;
probFeatUndet1 = numFeatUndet1 / numFeatTot;
probFeatUndet2 = numFeatUndet2 / numFeatTot;

%% merges/splits vs. features

%initialize lists of merges/splits
listOfMerges = [];
listOfSplits = [];

%go over all compound tracks
for iTrack = 1 : numTracks

    %get the sequence of events of this compound track
    seqOfEvents = tracks(iTrack).seqOfEvents;

    %find where this track has merges/splits
    indxMS = find(~isnan(seqOfEvents(:,4)));

    %go over these merges/splits
    for iMS = indxMS'

        %determine whether it's a merge (2) or a split (1)
        msType = seqOfEvents(iMS,2);

        %get the indices of the participating segments within the compound
        %track
        segmentsMS = seqOfEvents(iMS,3:4);

        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;

        %add this merge/split to the list of merges/splits
        if msType == 2 %merge
            listOfMerges = [listOfMerges; segmentsMS];
        else %split
            listOfSplits = [listOfSplits; segmentsMS];
        end

    end

end

%get total number of merges and splits
numMergesTot = size(listOfMerges,1);
numSplitsTot = size(listOfSplits,1);

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate average number of merges/splits per feature - this is the
%probability of a feature merging/splitting
probFeatMerge = aveMergePerFrame / aveFeatPerFrame;
probFeatSplit = aveSplitPerFrame / aveFeatPerFrame;

%% merges/splits vs. type

%get the types of segments participating in merges and splits
listOfMergeTypes = [trackSegmentType(listOfMerges(:,1)) ...
    trackSegmentType(listOfMerges(:,2))];
listOfSplitTypes = [trackSegmentType(listOfSplits(:,1)) ...
    trackSegmentType(listOfSplits(:,2))];

%sort the lists so that the "larger" (more dynamic) type is in the first
%column
listOfMergeTypes = sort(listOfMergeTypes,2,'descend');
listOfSplitTypes = sort(listOfSplitTypes,2,'descend');

%get number of merges/splits per type

%first approach:
%classify a merge/split type based on the more dynamic of its two segments 
%in terms of dynamics, linear > Brownian > confined Browian > undetermined 
%& length >= 5 > undetermined &length < 5
%for example, if a linear track segment merges with a Brownian track
%segment, the merge is classified as linear
numMergesLin1    = length(find(listOfMergeTypes(:,1)==3));
numMergesBrown1  = length(find(listOfMergeTypes(:,1)==2));
numMergesConf1   = length(find(listOfMergeTypes(:,1)==1));
numMergesUndet11 = length(find(listOfMergeTypes(:,1)==0));
numMergesUndet21 = length(find(listOfMergeTypes(:,1)==-1));
numSplitsLin1    = length(find(listOfSplitTypes(:,1)==3));
numSplitsBrown1  = length(find(listOfSplitTypes(:,1)==2));
numSplitsConf1   = length(find(listOfSplitTypes(:,1)==1));
numSplitsUndet11 = length(find(listOfSplitTypes(:,1)==0));
numSplitsUndet21 = length(find(listOfSplitTypes(:,1)==-1));

%second approach:
%get number of merges/splits involving a certain motion type regardless of
%the other segment type
%no ordering of dynamics in this case
indxMergesLin2    = find(any(listOfMergeTypes==3,2));
numMergesLin2     = length(indxMergesLin2);
indxMergesBrown2  = find(any(listOfMergeTypes==2,2));
numMergesBrown2   = length(indxMergesBrown2);
indxMergesConf2   = find(any(listOfMergeTypes==1,2));
numMergesConf2    = length(indxMergesConf2);
indxMergesUndet12 = find(any(listOfMergeTypes==0,2));
numMergesUndet12  = length(indxMergesUndet12);
indxMergesUndet22 = find(any(listOfMergeTypes==-1,2));
numMergesUndet22  = length(indxMergesUndet22);
indxSplitsLin2    = find(any(listOfSplitTypes==3,2));
numSplitsLin2     = length(indxSplitsLin2);
indxSplitsBrown2  = find(any(listOfSplitTypes==2,2));
numSplitsBrown2   = length(indxSplitsBrown2);
indxSplitsConf2   = find(any(listOfSplitTypes==1,2));
numSplitsConf2    = length(indxSplitsConf2);
indxSplitsUndet12 = find(any(listOfSplitTypes==0,2));
numSplitsUndet12  = length(indxSplitsUndet12);
indxSplitsUndet22 = find(any(listOfSplitTypes==-1,2));
numSplitsUndet22  = length(indxSplitsUndet22);

% % % %third approach:
% % % %classify merges/splits based on the types of the two track segments
% % % %involved
% % % numMergesLin3   = [length(find(listOfMergeTypes(indxMergesLin2,2)==3)) ...
% % %     length(find(listOfMergeTypes(indxMergesLin2,2)==2)) ...
% % %     length(find(listOfMergeTypes(indxMergesLin2,2)==1)) ...
% % %     length(find(listOfMergeTypes(indxMergesLin2,2)==0))];
% % % numMergesBrown3 = [length(find(listOfMergeTypes(indxMergesBrown2,1)==3)) ...
% % %     length(find(all(listOfMergeTypes(indxMergesBrown2,:)==2,2))) ...
% % %     length(find(listOfMergeTypes(indxMergesBrown2,2)==1)) ...
% % %     length(find(listOfMergeTypes(indxMergesBrown2,2)==0))];
% % % numMergesConf3  = [length(find(listOfMergeTypes(indxMergesConf2,1)==3)) ...
% % %     length(find(listOfMergeTypes(indxMergesConf2,1)==2)) ...
% % %     length(find(all(listOfMergeTypes(indxMergesConf2,:)==1,2))) ...
% % %     length(find(listOfMergeTypes(indxMergesConf2,2)==0))];
% % % numMergesUndet3 = [length(find(listOfMergeTypes(indxMergesUndet2,1)==3)) ...
% % %     length(find(listOfMergeTypes(indxMergesUndet2,1)==2)) ...
% % %     length(find(listOfMergeTypes(indxMergesUndet2,1)==1)) ...
% % %     length(find(listOfMergeTypes(indxMergesUndet2,1)==0))];
% % % numSplitsLin3   = [length(find(listOfSplitTypes(indxSplitsLin2,2)==3)) ...
% % %     length(find(listOfSplitTypes(indxSplitsLin2,2)==2)) ...
% % %     length(find(listOfSplitTypes(indxSplitsLin2,2)==1)) ...
% % %     length(find(listOfSplitTypes(indxSplitsLin2,2)==0))];
% % % numSplitsBrown3 = [length(find(listOfSplitTypes(indxSplitsBrown2,1)==3)) ...
% % %     length(find(all(listOfSplitTypes(indxSplitsBrown2,:)==2,2))) ...
% % %     length(find(listOfSplitTypes(indxSplitsBrown2,2)==1)) ...
% % %     length(find(listOfSplitTypes(indxSplitsBrown2,2)==0))];
% % % numSplitsConf3  = [length(find(listOfSplitTypes(indxSplitsConf2,1)==3)) ...
% % %     length(find(listOfSplitTypes(indxSplitsConf2,1)==2)) ...
% % %     length(find(all(listOfSplitTypes(indxSplitsConf2,:)==1,2))) ...
% % %     length(find(listOfSplitTypes(indxSplitsConf2,2)==0))];
% % % numSplitsUndet3 = [length(find(listOfSplitTypes(indxSplitsUndet2,1)==3)) ...
% % %     length(find(listOfSplitTypes(indxSplitsUndet2,1)==2)) ...
% % %     length(find(listOfSplitTypes(indxSplitsUndet2,1)==1)) ...
% % %     length(find(listOfSplitTypes(indxSplitsUndet2,1)==0))];

%get the fraction of merge/split types - this is the conditional
%probability of having a certain motion type if merging/splitting

%first approach 
probLinIfMerge1   = numMergesLin1   / numMergesTot;
probBrownIfMerge1 = numMergesBrown1 / numMergesTot;
probConfIfMerge1  = numMergesConf1  / numMergesTot;
probUndet1IfMerge1 = numMergesUndet11 / numMergesTot;
probUndet2IfMerge1 = numMergesUndet21 / numMergesTot;
probLinIfSplit1   = numSplitsLin1   / numSplitsTot;
probBrownIfSplit1 = numSplitsBrown1 / numSplitsTot;
probConfIfSplit1  = numSplitsConf1  / numSplitsTot;
probUndet1IfSplit1 = numSplitsUndet11 / numSplitsTot;
probUndet2IfSplit1 = numSplitsUndet21 / numSplitsTot;

%second approach 
probLinIfMerge2   = numMergesLin2   / numMergesTot;
probBrownIfMerge2 = numMergesBrown2 / numMergesTot;
probConfIfMerge2  = numMergesConf2  / numMergesTot;
probUndet1IfMerge2 = numMergesUndet12 / numMergesTot;
probUndet2IfMerge2 = numMergesUndet22 / numMergesTot;
probLinIfSplit2   = numSplitsLin2   / numSplitsTot;
probBrownIfSplit2 = numSplitsBrown2 / numSplitsTot;
probConfIfSplit2  = numSplitsConf2  / numSplitsTot;
probUndet1IfSplit2 = numSplitsUndet12 / numSplitsTot;
probUndet2IfSplit2 = numSplitsUndet22 / numSplitsTot;

% % % %third approach
% % % probLinIfMerge3   = numMergesLin3   / numMergesTot;
% % % probBrownIfMerge3 = numMergesBrown3 / numMergesTot;
% % % probConfIfMerge3  = numMergesConf3  / numMergesTot;
% % % probUndetIfMerge3 = numMergesUndet3 / numMergesTot;
% % % probLinIfSplit3   = numSplitsLin3   / numSplitsTot;
% % % probBrownIfSplit3 = numSplitsBrown3 / numSplitsTot;
% % % probConfIfSplit3  = numSplitsConf3  / numSplitsTot;
% % % probUndetIfSplit3 = numSplitsUndet3 / numSplitsTot;

%calculate the conditional probability of a feature undergoing a
%merge/split IF undergoing a certain motion type

%first approach
probMergeIfLin1   = probLinIfMerge1   * probFeatMerge / probFeatLin;
probMergeIfBrown1 = probBrownIfMerge1 * probFeatMerge / probFeatBrown;
probMergeIfConf1  = probConfIfMerge1  * probFeatMerge / probFeatConf;
probMergeIfUndet11 = probUndet1IfMerge1 * probFeatMerge / probFeatUndet1;
probMergeIfUndet21 = probUndet2IfMerge1 * probFeatMerge / probFeatUndet2;
probSplitIfLin1   = probLinIfSplit1   * probFeatSplit / probFeatLin;
probSplitIfBrown1 = probBrownIfSplit1 * probFeatSplit / probFeatBrown;
probSplitIfConf1  = probConfIfSplit1  * probFeatSplit / probFeatConf;
probSplitIfUndet11 = probUndet1IfSplit1 * probFeatSplit / probFeatUndet1;
probSplitIfUndet21 = probUndet2IfSplit1 * probFeatSplit / probFeatUndet2;

%second approach
probMergeIfLin2   = probLinIfMerge2   * probFeatMerge / probFeatLin;
probMergeIfBrown2 = probBrownIfMerge2 * probFeatMerge / probFeatBrown;
probMergeIfConf2  = probConfIfMerge2  * probFeatMerge / probFeatConf;
probMergeIfUndet12 = probUndet1IfMerge2 * probFeatMerge / probFeatUndet1;
probMergeIfUndet22 = probUndet2IfMerge2 * probFeatMerge / probFeatUndet2;
probSplitIfLin2   = probLinIfSplit2   * probFeatSplit / probFeatLin;
probSplitIfBrown2 = probBrownIfSplit2 * probFeatSplit / probFeatBrown;
probSplitIfConf2  = probConfIfSplit2  * probFeatSplit / probFeatConf;
probSplitIfUndet12 = probUndet1IfSplit2 * probFeatSplit / probFeatUndet1;
probSplitIfUndet22 = probUndet2IfSplit2 * probFeatSplit / probFeatUndet2;

% % % %calculate the conditional probability of 2 features undergoing a
% % % %merge/split IF each is undergoing a certain motion type
% % % 
% % % %third approach
% % % probFeatAll = [probFeatLin probFeatBrown probFeatConf probFeatUndet];
% % % probMergeIfLin3   = probLinIfMerge3   * probFeatMerge / probFeatLin   ./ probFeatAll;
% % % probMergeIfBrown3 = probBrownIfMerge3 * probFeatMerge / probFeatBrown ./ probFeatAll;
% % % probMergeIfConf3  = probConfIfMerge3  * probFeatMerge / probFeatConf  ./ probFeatAll;
% % % probMergeIfUndet3 = probUndetIfMerge3 * probFeatMerge / probFeatUndet ./ probFeatAll;
% % % probSplitIfLin3   = probLinIfSplit3   * probFeatSplit / probFeatLin   ./ probFeatAll;
% % % probSplitIfBrown3 = probBrownIfSplit3 * probFeatSplit / probFeatBrown ./ probFeatAll;
% % % probSplitIfConf3  = probConfIfSplit3  * probFeatSplit / probFeatConf  ./ probFeatAll;
% % % probSplitIfUndet3 = probUndetIfSplit3 * probFeatSplit / probFeatUndet ./ probFeatAll;

%% time from merges to splits and vice versa

% for iType = 0 : 1
%
%     if iType == 1
%         trackType = 'Lin';
%     else
%         trackType =  'Brown';
%     end
%
%     %initialize some variables
%     numTrackNoMS = 0;
%     numTrackOnlyM = 0;
%     numTrackOnlyS = 0;
%     numTrackMS = 0;
%     timeMerge2Split = [];
%     timeSplit2Merge = [];
%     eval(['indxTracks = indx' trackType ';']);
%
%     %go over all tracks of this type ...
%     for iTrack = indxTracks'
%
%         %get track's sequence of events
%         seqOfEvents = tracks(iTrack).seqOfEvents;
%
%         %find merge and split times in track
%         msTime = seqOfEvents(~isnan(seqOfEvents(:,4)),[1 2]);
%         msTime(msTime(:,2)==1,1) = -msTime(msTime(:,2)==1,1);
%         msTime = msTime(:,1);
%
%         %take action based on whether there are merges and splits
%         if isempty(msTime) %if there are no merges and no splits
%
%             %add one to counter of tracks without merges and splits
%             numTrackNoMS = numTrackNoMS + 1;
%
%         elseif all(msTime > 0) %if there are merges but no splits
%
%             %add one to counter of tracks with only merges
%             numTrackOnlyM = numTrackOnlyM + 1;
%
%         elseif all(msTime < 0) %if there are splits but no merges
%
%             %add one to counter of tracks with only splits
%             numTrackOnlyS = numTrackOnlyS + 1;
%
%         else %if there are both merges and splits
%
%             %add one to counter of tracks with both merges and
%             %splits
%             numTrackMS = numTrackMS + 1;
%
%             %get number of merges and splits to consider
%             numMS = length(msTime);
%
%             indxMS = 1;
%             while indxMS < numMS
%
%                 %get index of initial event and its type
%                 typeEvent1 = sign(msTime(indxMS));
%                 indxEvent1 = indxMS;
%
%                 %keep increasing indxMS until you find an event different
%                 %from initial event
%                 typeEvent2 = typeEvent1;
%                 while typeEvent2 == typeEvent1 && indxMS < numMS
%                     indxMS = indxMS + 1;
%                     typeEvent2 = sign(msTime(indxMS));
%                 end
%                 indxEvent2 = indxMS;
%
%                 %keep increasing indxMS until you find an event again
%                 %similar to initial event
%                 typeEvent3 = typeEvent2;
%                 while typeEvent3 ~= typeEvent1 && indxMS < numMS
%                     indxMS = indxMS + 1;
%                     typeEvent3 = sign(msTime(indxMS));
%                 end
%                 if typeEvent3 ~= typeEvent2
%                     indxMS = indxMS - 1;
%                 end
%                 indxEvent3 = indxMS;
%
%                 if typeEvent2 ~= typeEvent1
%
%                     %calculate total number of combinations for calculating
%                     %time between merges and splits (or vice versa)
%                     numCombination = (indxEvent2-indxEvent1)*(indxEvent3-indxEvent2+1);
%
%                     %go over all combinations and calculate time
%                     indxComb = 0;
%                     timeBetweenMS = zeros(numCombination,1);
%                     for indx1 = 1 : indxEvent2-indxEvent1
%                         for indx2 = 1 : indxEvent3-indxEvent2+1
%                             indxComb = indxComb + 1;
%                             timeBetweenMS(indxComb) = abs(msTime(indx2+indxEvent2-1)) - ...
%                                 abs(msTime(indx1+indxEvent1-1));
%                         end
%                     end
%
%                     %store the times and their weights based on whether we
%                     %looked at a merge to split or a split to merge
%                     if typeEvent1 > 0
%                         timeMerge2Split = [timeMerge2Split; [timeBetweenMS ...
%                             (1/numCombination)*ones(numCombination,1)]];
%                     else
%                         timeSplit2Merge = [timeSplit2Merge; [timeBetweenMS ...
%                             (1/numCombination)*ones(numCombination,1)]];
%                     end
%
%                 end
%
%                 %update indxMS to look at next merge-to-split or
%                 %split-to-merge event
%                 indxMS = indxEvent2;
%
%             end %(while indxMS <= numMS - 1)
%
%         end %(if isempty(msTime) ... elseif ...)
%
%     end %(for iTrack = 1 : indxTracks')
%
%     %store track numbers and distributions
%     eval(['numTracks' trackType ' = [numTrackMS numTrackOnlyM numTrackOnlyS numTrackNoMS];'])
%     eval(['timeMerge2Split' trackType ' = timeMerge2Split;'])
%     eval(['timeSplit2Merge' trackType ' = timeSplit2Merge;'])
%
% end
    
%% output

%statistics per category of motion
statsPerCat = [fracSegmentsLin probFeatLin probMergeIfLin1 probSplitIfLin1 ...
    probMergeIfLin2 probSplitIfLin2; ...
    fracSegmentsBrown probFeatBrown probMergeIfBrown1 probSplitIfBrown1 ...
    probMergeIfBrown2 probSplitIfBrown2; ...
    fracSegmentsConf probFeatConf probMergeIfConf1 probSplitIfConf1 ...
    probMergeIfConf2 probSplitIfConf2; ...
    fracSegmentsUndet1 probFeatUndet1 probMergeIfUndet11 probSplitIfUndet11 ...
    probMergeIfUndet12 probSplitIfUndet12; ...
    fracSegmentsUndet2 probFeatUndet2 probMergeIfUndet21 probSplitIfUndet21 ...
    probMergeIfUndet22 probSplitIfUndet22];
    
% % % statsPerCat = [fracSegmentsLin probFeatLin probMergeIfLin1 probSplitIfLin1 ...
% % %     probMergeIfLin2 probSplitIfLin2 probMergeIfLin3 probSplitIfLin3; ...
% % %     fracSegmentsBrown probFeatBrown probMergeIfBrown1 probSplitIfBrown1 ...
% % %     probMergeIfBrown2 probSplitIfBrown2 probMergeIfBrown3 probSplitIfBrown3; ...
% % %     fracSegmentsConf probFeatConf probMergeIfConf1 probSplitIfConf1 ...
% % %     probMergeIfConf2 probSplitIfConf2 probMergeIfConf3 probSplitIfConf3; ...
% % %     fracSegmentsUndet probFeatUndet probMergeIfUndet1 probSplitIfUndet1 ...
% % %     probMergeIfUndet2 probSplitIfUndet2 probMergeIfUndet3 probSplitIfUndet3];

%additional general statistics
statsGeneral = [aveFeatPerFrame numTracks numTrackSegments probFeatMerge probFeatSplit];

%number of events per category, to get an idea of how many observations
%there are
numEventsPerCat = [numMergesLin1 numSplitsLin1 numMergesLin2 numSplitsLin2; ...
    numMergesBrown1 numSplitsBrown1 numMergesBrown2 numSplitsBrown2; ...
    numMergesConf1 numSplitsConf1 numMergesConf2 numSplitsConf2; ...
    numMergesUndet11 numSplitsUndet11 numMergesUndet12 numSplitsUndet12; ...
    numMergesUndet21 numSplitsUndet21 numMergesUndet22 numSplitsUndet22];

% % % numEventsPerCat = [numMergesLin1 numSplitsLin1 numMergesLin2 numSplitsLin2 ...
% % %     numMergesLin3 numSplitsLin3; ...
% % %     numMergesBrown1 numSplitsBrown1 numMergesBrown2 numSplitsBrown2 ...
% % %     numMergesBrown3 numSplitsBrown3; ...
% % %     numMergesConf1 numSplitsConf1 numMergesConf2 numSplitsConf2 ...
% % %     numMergesConf3 numSplitsConf3; ...
% % %     numMergesUndet1 numSplitsUndet1 numMergesUndet2 numSplitsUndet2 ...
% % %     numMergesUndet3 numSplitsUndet3];

% msTimeInfo.brown.numTracks = numTracksBrown;
% msTimeInfo.brown.timeMerge2Split = timeMerge2SplitBrown;
% msTimeInfo.brown.timeSplit2Merge = timeSplit2MergeBrown;
% msTimeInfo.linear.numTracks = numTracksLin;
% msTimeInfo.linear.timeMerge2Split = timeMerge2SplitLin;
% msTimeInfo.linear.timeSplit2Merge = timeSplit2MergeLin;
msTimeInfo = [];

%% ~~~ the end ~~~


%% trial stuff

% % % for iTrack = 1 : numTracks
% % %     trackSegmentClassTmp = diffAnalysisRes(iTrack).classification;
% % %     if any(trackSegmentClassTmp(:,1) == 1)
% % %         trackClass = 3;
% % %     elseif any(trackSegmentClassTmp(:,2) == 2)
% % %         trackClass = 2;
% % %     elseif any(trackSegmentClassTmp(:,2) == 1)
% % %         trackClass = 1;
% % %     else
% % %         trackClass = 0;
% % %     end
% % %     trackSegmentType(trackStartRow(iTrack):trackStartRow(iTrack)+...
% % %         numSegments(iTrack)-1,1) = trackClass;
% % % end
% % % indxLin = find(trackSegmentType==3);
% % % indxBrown = find(trackSegmentType==2);
% % % indxConf = find(trackSegmentType==1);
% % % indxUndet = find(trackSegmentType==0);


