function plotCompTrackAmp(trackedFeatureInfo)
%PLOTCOMPTRACKAMP plots the intensities along a compound track, indicating merges, splits and gaps
%
%SYNOPSIS plotCompTrackAmp(trackedFeatureInfo)
%
%INPUT  trackedFeatureInfo: Output of trackCloseGapsKalman for one track:
%                           Contains the fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%
%OUTPUT The plot.
%
%REMARKS gaps are dotted black lines, splits are dash-dotted black lines
%and merges are dashed lines
%
%Khuloud Jaqaman, May 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--plotTracks2D: Incorrect number of input arguments!');
    return
end

%extract information from input
seqOfEvents = trackedFeatureInfo.seqOfEvents;
tracksCoordAmpCG = trackedFeatureInfo.tracksCoordAmpCG;

%get first frame, last frame and number of frames
firstFrame = seqOfEvents(1,1);
lastFrame = seqOfEvents(end,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract amplitudes from input
ampSequence = tracksCoordAmpCG(:,4:8:end);

%find number of segments making compound track
numSegments = size(ampSequence,1);

%make new figure and hold on to it
figure, hold on

%plot amplitudes as dotted black lines, closing gaps
for i = 1 : numSegments
    indx = find(~isnan(ampSequence(i,:)));
    plot(indx+firstFrame-1,ampSequence(i,indx),'k:');
end

%plot amplitudes in color leaving gaps as blank (so that they appear as
%dotted lines in the final figure)
plot((firstFrame:lastFrame)',ampSequence','marker','.');

%find merges and splits
indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';

%go over all splits
for iSplit = indxSplit

    %get time of splitting
    timeSplit = seqOfEvents(iSplit,1);

    %determine index of starting track
    rowS = seqOfEvents(iSplit,3);

    %determine index of splitting track
    rowSp = seqOfEvents(iSplit,4);

    %plot split as a black dash-dotted line
    plot([timeSplit-1 timeSplit],[ampSequence(rowSp,timeSplit-1) ...
        ampSequence(rowS,timeSplit)],'k-.') 
    
end

%go over all merges
for iMerge = indxMerge

    %get time of merging
    timeMerge = seqOfEvents(iMerge,1);

    %determine index of ending track
    rowE = seqOfEvents(iMerge,3);

    %determine index of merging track
    rowM = seqOfEvents(iMerge,4);

    %plot merge as a black dashed line
    plot([timeMerge-1 timeMerge],[ampSequence(rowE,timeMerge-1) ...
        ampSequence(rowM,timeMerge)],'k--') 
    
end

%hold off of figure
hold off


%%%%% ~~ the end ~~ %%%%%

