function [costVec,errFlag] = costVecLinkMS(trackedFeatMS,timeSplit,...
    timeMerge,costMatParams,gapCloseParam)
%COSTVECLINKMS provides the costs of linking merging tracks to a splitting track
%
%SYNOPSIS [costVec,errFlag] = costVecLinkMS(trackedFeatMS,timeSplit,...
%    timeMerge,costMatParams,gapCloseParam)
%
%INPUT  trackedFeatMS: The positions and amplitudes of the spltting track
%                      (in 1st row) and the possible merging tracks (in the
%                      rest of the rows). Each row consists of 
%                      [x1 y1 a1 dx1 dy1 da1 x2 y2 a2 dx2 dy2 da2 ...]
%                      in image coordinate system (coordinates in
%                      pixels). NaN is used to indicate time points 
%                        where the track does not exist.
%       timeSplit    : Time of splitting.
%       timeMerge    : Time(s) of merging.
%       costMatParams: Structure with the fields:
%             .trackStats : Structure with the following fields:
%                   .dispSqLambda: Parameters of the exponential distribution
%                                  that describes the displacement of a feature
%                                  between two time points.
%                   .ampDiffStd  : Standard deviations of the change in a feature's 
%                                  amplitude between two time points.
%             .cutoffProbD: Cumulative probability of a square diplacement
%                           beyond which linking is not allowed.
%             .cutoffProbA: Cumulative probability of an amplitude
%                           difference beyond which linking is not allowed.
%       gapCloseParam: Structure containing variables needed for gap closing.
%                      Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to it.
%
%OUTPUT costVec      : Vector of costs.
%       errFlag      : 0 if function executes normally, 1 otherwise.
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
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costVec = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costVecLinkMS')
    disp('--costVecLinkMS: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

dispSqLambda = costMatParams.trackStats.dispSqLambda;
ampDiffStd = costMatParams.trackStats.ampDiffStd;
cutoffProbD = costMatParams.cutoffProbD;
cutoffProbA = costMatParams.cutoffProbA;
timeWindow = gapCloseParam.timeWindow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost vector calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the maximum squared displacement that allows linking 2 features
maxDispSq = expinv(cutoffProbD,dispSqLambda);

%find the maximum squared amplitude change that allows linking 2 features
maxAmpDiffSq = (norminv(cutoffProbA,0,ampDiffStd)).^2;

%calculate the additive constant for each cost as a function of time gap
addConst = log(ampDiffStd) - log(dispSqLambda);

%get the number of possible merging tracks
numPrevMerge = size(trackedFeatMS,1) - 1;

%get the x,y-coordinates and amplitude at time of splitting
xCoordStart = trackedFeatMS(1,(timeSplit-1)*6+1);
yCoordStart = trackedFeatMS(1,(timeSplit-1)*6+2);
ampStart    = trackedFeatMS(1,(timeSplit-1)*6+3);

%get the x,y-coordinates and amplitudes before merging
xCoordEnd = zeros(numPrevMerge,1);
yCoordEnd = zeros(numPrevMerge,1);
ampEnd    = zeros(numPrevMerge,1);
for i=1:numPrevMerge
    xCoordEnd(i) = trackedFeatMS(i+1,(timeMerge(i)-2)*6+1);
    yCoordEnd(i) = trackedFeatMS(i+1,(timeMerge(i)-2)*6+2);
    ampEnd(i)    = trackedFeatMS(i+1,(timeMerge(i)-2)*6+3);
end

%initialize cost vector
costVec  = Inf*ones(numPrevMerge,1);

for i=1:numPrevMerge %go over all possible merges

    %subtract splitting time from time just before merging
    timeGap = timeSplit - timeMerge(i) + 1;

    %if the time gap between merging and splitting is within the acceptable limits ...
    if timeGap > 1 && timeGap <= timeWindow

        %calculate the square distance between the end
        %point and starting point
        dispSq = (xCoordEnd(i)-xCoordStart)^2 + ...
            (yCoordEnd(i)-yCoordStart)^2;

        %calculate the square difference between the amplitude at
        %the beginning of track j and the end of track i
        ampDiffSq = (ampEnd(i) - ampStart)^2;

        %if this is a possible link ...
        if dispSq < maxDispSq(timeGap) && ...
                ampDiffSq < maxAmpDiffSq(timeGap)

            %assign the cost for this pair
            costVec(i) = dispSqLambda(timeGap)*dispSq + ...
                ampDiffSq/2/ampDiffStd(timeGap)^2 + ...
                addConst(timeGap);

        end %(if dispSq < maxDispSq(timeGap) && ampDIffSq < maxAmpDiffSq(timeGap))

    end %(if timeGap > 1 || timeGap <= timeWindow)

end %(for i=1:numPrevMerge)


%%%%% ~~ the end ~~ %%%%%

