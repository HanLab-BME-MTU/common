function [costMat,noLinkCost,trackStartTime,trackEndTime,indxMerge,numMerge,...
    indxSplit,numSplit,errFlag] = costMatCloseGaps(trackedFeatInfo,...
    trackStartTime,indxStart,trackEndTime,indxEnd,costMatParams,gapCloseParam)
%COSTMATCLOSEGAPS provides a cost matrix for closing gaps (including merging and splitting) based on tracked feature statistics
%
%SYNOPSIS [costMat,noLinkCost,trackStartTime,trackEndTime,indxMerge,numMerge,...
%    indxSplit,numSplit,errFlag] = costMatCloseGaps(trackedFeatInfo,...
%    trackStartTime,indxStart,trackEndTime,indxEnd,costMatParams,gapCloseParam)
%
%INPUT  trackedFeatInfo: The positions and amplitudes of the tracked
%                        features from linkFeaturesTp2Tp. 
%                        Number of rows = number of tracks.
%                        Number of columns = 6*number of time points. 
%                        Each row consists of 
%                        [x1 y1 a1 dx1 dy1 da1 x2 y2 a2 dx2 dy2 da2 ...]
%                        in image coordinate system (coordinates in
%                        pixels). NaN is used to indicate time points 
%                        where the track does not exist.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       costMatParams  : Structure with the fields:
%             .trackStats : Structure with the following fields:
%                   .dispSqLambda: Parameters of the exponential distribution
%                                  that describes the displacement of a feature
%                                  between two time points.
%                   .ampDiffStd  : Standard deviations of the change in a feature's 
%                                  amplitude between two time points.
%             .cutCProb1  : Cumulative probability of a square diplacement
%                           or amplitude difference beyond which linking
%                           between an end and a start is not allowed.
%             .cutCProb2  : Cumulative probability of a square displacement
%                           or amplitude difference beyond which merging
%                           and splitting are not allowed.
%             .noLnkPrctl : Percentile used to calculate the cost of
%                           linking a feature to nothing. Use -1 if you do
%                           not want to calculate this cost.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to it.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%
%OUTPUT costMat       : Cost matrix.
%       noLinkCost    : Cost of linking a feature to nothing, as derived
%                       from the distribution of costs.
%       trackStartTime: Starting time of tracks whose starts are considered
%                       for gap closing.
%       trackEndTime  : Ending time of tracks whose ends are considered for
%                       gap closing.
%       indxMerge     : Index of tracks that have possibly merged with
%                       tracks that end before the last time points.
%       numMerge      : Number of such tracks.
%       indxSplit     : Index of tracks from which tracks that begin after
%                       the first time point might have split.
%       numSplit      : Number of such tracks.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%REMARKS the costs are given by ...
%
%The cost for linking the end of track i at time point t to the 
%beginning of track j at time point t' (t' > t) is given by
%-log[p(dI)p(dispSq)], where dI = I(j;t') - I(i,t) is assumed to be 
%normally distributed with mean 0 and standard deviation ampDiffStd(t'-t) 
%(supplied by user), and dispSq = square of distance between end position 
%of i at t and beginning position of j at t' is assumed to be exponentially
%distributed with parameter dispSqLambda(t'-t) (supplied by user).
%
%The cost for merging and splitting is given by -log[p(dI)p(dispSq)], where
%dispSq and p(dispSq) are as described above and, in the case of merging,
%dI = I(feature after merging) - sum(I(features before merging)), and, 
%in the case of splitting, 
%dI = sum(I(features after splitting)) - I(feature before splitting). 
%dI is thus normally distributed with mean zeros and standard deviation
%sqrt(2)*ampDiffStd(1).
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costMat = [];
noLinkCost = [];
indxMerge = [];
numMerge = [];
indxSplit = [];
numSplit = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatCloseGaps')
    disp('--costMatCloseGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

timeWindow = gapCloseParam.timeWindow;
mergeSplit = gapCloseParam.mergeSplit;
dispSqLambda = costMatParams.trackStats.dispSqLambda;
ampDiffStd = costMatParams.trackStats.ampDiffStd;
cutCProb1 = costMatParams.cutCProb1;
if mergeSplit
    cutCProb2 = costMatParams.cutCProb2;
end
noLnkPrctl = costMatParams.noLnkPrctl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the maximum squared displacement that allows linking an end to a start
maxDispSq = expinv(cutCProb1,dispSqLambda);

%find the maximum squared amplitude change that allows linking an end to a start
maxAmpDiffSq = (norminv(cutCProb1,0,ampDiffStd)).^2;

%calculate the additive constant for each cost as a function of time gap
addConst = log(ampDiffStd) - log(dispSqLambda);

%get number of tracks formed by initial linking
numTracks = size(trackedFeatInfo,1);

%get the number of tracks whose starts are to be considered
m = length(indxStart);

%get the number of tracks whose ends are to be considered
n = length(indxEnd);

%get the x,y-coordinates and amplitudes of features at the starts of
%their tracks
xCoordStart = zeros(numTracks,1);
yCoordStart = zeros(numTracks,1);
ampStart = zeros(numTracks,1);
for i=1:numTracks
    xCoordStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*6+1);
    yCoordStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*6+2);
    ampStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*6+3);
end

%get the x,y-coordinates and amplitudes of features at the ends of
%their tracks
xCoordEnd = zeros(numTracks,1);
yCoordEnd = zeros(numTracks,1);
ampEnd = zeros(numTracks,1);
for i=1:numTracks
    xCoordEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*6+1);
    yCoordEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*6+2);
    ampEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*6+3);
end

%remove tracks that start at the first time point and get the total number of
%tracks whose starts are to be considered
trackStartTime = trackStartTime(indxStart);
xCoordStart = xCoordStart(indxStart);
yCoordStart = yCoordStart(indxStart);
ampStart = ampStart(indxStart);

%remove tracks that end at the last time point and get the total number of
%tracks whose ends are to be considered
trackEndTime = trackEndTime(indxEnd);
xCoordEnd = xCoordEnd(indxEnd);
yCoordEnd = yCoordEnd(indxEnd);
ampEnd = ampEnd(indxEnd);

indx1 = []; %row number in cost matrix
indx2 = []; %column number in cost matrix
cost  = []; %cost value

%costs for closing gaps due to features going out-of-focus

for j=1:m %go over all starts
    for i=1:n %go over all ends

        %subtract starting time from ending time
        timeGap = trackStartTime(j) - trackEndTime(i);

        %if the time gap between the end of track i and the start of track
        %j is within the acceptable limits ...
        if timeGap > 1 && timeGap <= timeWindow

            %calculate the square distance between the end
            %point and starting point
            dispSq = (xCoordEnd(i)-xCoordStart(j))^2 + ...
                (yCoordEnd(i)-yCoordStart(j))^2;

            %calculate the square difference between the amplitude at
            %the beginning of track j and the end of track i
            ampDiffSq = (ampEnd(i) - ampStart(j))^2;

            %if this is a possible link ...
            if dispSq < maxDispSq(timeGap) && ...
                    ampDiffSq < maxAmpDiffSq(timeGap)

                %assign the cost for this pair of tracks
                indx1 = [indx1; i]; %row number
                indx2 = [indx2; j]; %column number
                cost = [cost; dispSqLambda(timeGap)*dispSq + ...
                    ampDiffSq/2/ampDiffStd(timeGap)^2 + ...
                    addConst(timeGap)]; %cost

            end %(if dispSq < maxDispSq(timeGap) && ampDIffSq < maxAmpDiffSq(timeGap))

        end %(if timeGap > 1 || timeGap <= timeWindow)

    end %(for i=1:n)
end %(for j=1:m)

%determine the cost of birth and death
costBD = max(cost(:)) + 1;

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting

%if merging and splitting are to be considered ...
if mergeSplit

    %calculate additional additive constants to cost
    halfLog2 = log(2)/2;
    addConst2 = addConst(1) - log(dispSqLambda(1));
    
    %get the maximum squared displacement that allows merging and splitting
    maxDispSqMS = expinv(cutCProb2,dispSqLambda);

    %get maximum allowed intensity variation when merging or
    %splitting
    maxAmpDiffSqMS = (norminv(cutCProb2,0,1.4142*ampDiffStd)).^2;
    
    %costs of merging
    
    %go over all track ending times
    for endTime = min(trackEndTime):max(trackEndTime)

        %find tracks that end at this time point
        tracksToConsider = find(trackEndTime==endTime);
        numTracksToConsider = length(tracksToConsider);

        %first consider the case where the ending track merges with another
        %existing track. This happens when the difference in the intensities
        %of the two merging features is smaller than the change in intensity
        %from one time point to the next

        %get index indicating time of merging
        timeIndx  = endTime*6;

        for j=tracksToConsider' %go over all ends considered
            for i=1:numTracks %go over all tracks

                %get position and amplitude of merged feature
                xCoordMid = trackedFeatInfo(i,timeIndx+1);
                yCoordMid = trackedFeatInfo(i,timeIndx+2);
                ampMidT1 = trackedFeatInfo(i,timeIndx+3);

                %get amplitude of feature that merged with the
                %ending track, at time point before merging
                ampMidT0 = trackedFeatInfo(i,timeIndx-3);

                %calculate the square distance between the ending track
                %and the point of merging
                %dispSq = NaN if the track does not exist at the
                %time of merging (prohibiting a link)
                dispSq = (xCoordEnd(j)-xCoordMid)^2 ...
                    + (yCoordEnd(j)-yCoordMid)^2;

                %calculate the square difference between the amplitude
                %after merging and the sum of amplitudes before merging
                %ampDiffSq = NaN if the track does not exist at the
                %time of merging or the time point before it (prohibiting a
                %link).
                ampDiffSq = (ampEnd(j) + ampMidT0 - ampMidT1)^2;

                %if this is a possible link ...
                if dispSq < maxDispSqMS(1) && ...
                        ampDiffSq < maxAmpDiffSqMS(1)

                    %increase the "merge index" by one
                    numMerge = numMerge + 1;

                    %save the merging track's number
                    indxMerge = [indxMerge; i];

                    %calculate the cost of merging
                    indx1 = [indx1; j]; %row number
                    indx2 = [indx2; numMerge+m]; %column number

                    cost = [cost; dispSqLambda(1)*dispSq + ...
                        ampDiffSq/4/ampDiffStd(1)^2 + ...
                        addConst(1) + halfLog2]; %cost of merging

                    altCostMerge = [altCostMerge; 1 + halfLog2 + ... %cost of not mering
                        (ampMidT0-ampMidT1)^2/4/ampDiffStd(1)^2 + addConst(1)];

                end %(if dispSq < maxDispSqMS(1) && ampDiffSq < maxAmpDiffSqMS)

            end %(for i=1:numTracks)
        end %(for j=tracksToConsider')

%         %next consider the case where two ending tracks merge and form a 
%         %third track that starts right after they end. This happens when 
%         %the difference in the intensities of the two merging features 
%         %is larger than the change in intensity from one time point 
%         %to the next
% 
%         %find tracks that start right after these tracks end
%         startTracksToConsider = find(trackStartTime==endTime+1);
%         
%         for j1=1:numTracksToConsider %go over all pairs of ends considered
%             for k1=j1+1:numTracksToConsider
% 
%                 %get the end numbers
%                 j = tracksToConsider(j1);
%                 k = tracksToConsider(k1);
% 
%                 for i=startTracksToConsider'
% 
%                     %calculate the square distance between the end
%                     %point of j and starting point of i
%                     dispSq1 = (xCoordEnd(j)-xCoordStart(i))^2 + ...
%                         (yCoordEnd(j)-yCoordStart(i))^2;
% 
%                     %calculate the square distance between the end
%                     %point of k and starting point of i
%                     dispSq2 = (xCoordEnd(k)-xCoordStart(i))^2 + ...
%                         (yCoordEnd(k)-yCoordStart(i))^2;
% 
%                     %calculate the square difference between the amplitude
%                     %of i (after merging) and the sum of amplitudes of j
%                     %and k (before merging)
%                     ampDiffSq = (ampEnd(j) + ampEnd(k) - ampStart(i))^2;
% 
%                     %if this is a possible link ...
%                     if dispSq1 < maxDispSqMS(1) && ...
%                             dispSq2 < maxDispSqMS(1) && ...
%                             ampDiffSq < maxAmpDiffSqMS(1)
% 
%                         
% 
%                     end %(if dispSq1 < maxDispSqMS(1) && ...)
% 
%                 end %(for i=startTracksToConsider')
% 
%             end %(for k1=j1+1:numTracksToConsider)
%         end %(for j1=1:numTracksToConsider)

    end %(for endTime = min(trackEndTime):max(trackEndTime))

    %costs of splitting

    %go over all track starting times
    for startTime = min(trackStartTime):max(trackStartTime)
    
        %find tracks that start at this time point
        tracksToConsider = find(trackStartTime==startTime);

        %first consider the case where the starting track splits from another 
        %existing track. This happens when the difference in the intensities
        %of the two merged features is smaller than the change in intensity
        %from one time point to the next

        for j=tracksToConsider' %go over all starts considered
            for i=1:numTracks %go over all tracks

                %get indx indicating time of splitting
                timeIndx  = (trackStartTime(j)-2)*6;

                %get position and amplitude of feature before splitting
                xCoordMid = trackedFeatInfo(i,timeIndx+1);
                yCoordMid = trackedFeatInfo(i,timeIndx+2);
                ampMidT1 = trackedFeatInfo(i,timeIndx+3);

                %get amplitude of feature that split at time point
                %after splitting
                ampMidT0 = trackedFeatInfo(i,timeIndx+9);

                %calculate the square distance between the starting track
                %and the point of splitting
                dispSq = (xCoordStart(j)-xCoordMid)^2 ...
                    + (yCoordStart(j)-yCoordMid)^2;

                %calculate the square difference between the amplitude
                %after merging and the sum of amplitudes before merging
                ampDiffSq = (ampStart(j) + ampMidT0 - ampMidT1)^2;

                %if this is a possible link ...
                if dispSq < maxDispSqMS(1) && ...
                        ampDiffSq < maxAmpDiffSqMS(1)

                    %increase the "split index" by one
                    numSplit = numSplit + 1;

                    %save the splitting track's number
                    indxSplit = [indxSplit; i];

                    %calculate the cost of splitting
                    indx1 = [indx1; numSplit+n]; %row number
                    indx2 = [indx2; j]; %column number

                    cost = [cost; dispSqLambda(1)*dispSq + ...
                        ampDiffSq/4/ampDiffStd(1)^2 + ...
                        addConst(1) + halfLog2]; %cost of splitting

                    altCostSplit = [altCostSplit; 1 + halfLog2 + ... %cost of not mering
                        (ampMidT0-ampMidT1)^2/4/ampDiffStd(1)^2 + addConst(1)];

                end %(if dispSq < maxDispSqMS(1) && ampDiffSq < maxAmpDiffSqMS)

            end %(for i=1:numTracks)
        end %(for j=tracksToConsider')

    end %(for startTime = min(trackStartTime):max(trackStartTime))
        
end %(if mergeSplit)

%create temporary cost matrix
n1 = n+numSplit;
m1 = m+numMerge;
costMat = sparse(indx1,indx2,cost,n1,m1);

%append cost matrix to allow births and deaths ...

% %add block allowing deaths
% indx1 = [indx1; [1:n1]'];
% indx2 = [indx2; [m1+1:m1+n1]'];
% cost = [cost; costBD*ones(n,1); altCostSplit];
% 
% %add block allowing births
% indx1 = [indx1; [n1+1:n1+m1]'];
% indx2 = [indx2; [1:m1]'];
% cost = [cost; costBD*ones(m,1); altCostMerge];
% 
% %add last block (which has no practical use, but alas takes the most space...)
% indx1 = [indx1; repmat([n1+1:n1+m1]',n1,1)];
% indx2 = [indx2; reshape(repmat([m1+1:m1+n1],m1,1),n1*m1,1)];
% cost = [cost; ones(n1*m1,1)];
% 
% %create cost matrix that allows for births and deaths
% costMat = sparse(indx1,indx2,cost);

%determine noLinkCost
if noLnkPrctl ~= -1
    noLinkCost = prctile(cost,noLnkPrctl);
end


%%%%% ~~ the end ~~ %%%%%

