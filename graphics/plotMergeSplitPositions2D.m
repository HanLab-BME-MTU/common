function plotMergeSplitPositions2D(tracks,mergesInfo,splitsInfo,...
    mergesInfoSpace,splitsInfoSpace,figureName,image,mergesInfoSpaceExt,...
    splitsInfoSpaceExt)
%PLOTMERGESPLITPOSITIONS2D generates a spatial map of merge and split locations
%
%SYNPOSIS plotMergeSplitPositions2D(tracks,mergesInfo,splitsInfo,...
%    mergesInfoSpace,splitsInfoSpace,figureName,image,mergesInfoSpaceExt,...
%    splitsInfoSpaceExt)
%
%INPUT  tracks         : Output of trackCloseGapsKalman.
%       mergesInfo     : 2D array where first column indicates track number,
%                        second column indicates track type (1 linear, 0 o.w.),
%                        third column indicates number of merges, and 
%                        subsequent columns indicate merge times. Note
%                        that track type is not relevant, but is simply
%                        the output of findMergesSplits.
%       splitsInfo     : 2D array where first column indicates track number,
%                        second column indicates track type (1 linear, 0 o.w.),
%                        third column indicates number of splits, and 
%                        subsequence columns indicate split times. Note
%                        that track type is not relevant, but is simply
%                        the output of findMergesSplits.
%       mergesInfoSpace: 2D array that is a continuation of mergesInfoTime,
%                        storing the (x,y)-coordinates of each merge.
%                        Every row corresponds to the same row in
%                        mergesInfo. Every merge gets 2 columns, 1 for x
%                        and 1 for y.
%       splitsInfoSpace: 2D array that is a continuation of splitsInfoTime,
%                        storing the (x,y)-coordinates of each split.
%                        Every row corresponds to the same row in
%                        splitsInfo. Every split gets 2 columns, 1 for x
%                        and 1 for y.
%       figureName     : Figure name.
%                        Optional. Default: None.
%       image          : Image to overlay spatial map on.
%                        Optional. Default: [].
%       mergesInfoSpaceExt: 3D array that is a continuation of mergesInfoSpace
%                        but replicated in the 3rd dimension, storing the
%                        (x,y)-coordinates of the two particles before
%                        merging. First index of 3rd dimension is for
%                        continuing track segment, second index is for
%                        terminatating track segment.
%                        Optional. Default: []. If not supplied, only
%                        coordinate after merging (mergesInfoSpace) will be
%                        plotted.
%       splitsInfoSpaceExt: 3D array that is a continuation of splitsInfoSpace
%                        but replicated in the 3rd dimension, storing the
%                        (x,y)-coordinates of the two particles after
%                        merging. First index of 3rd dimension is for
%                        continuing track segment, second index is for
%                        initiating track segment.
%                        Optional. Default: []. If not supplied, only
%                        coordinate before splitting (splitsInfoSpace) will
%                        be plotted.
%
%REMARKS Merges are plotted in red, splits in black, while trajectories
%without merges or splits are plotted in gray.
%
%Khuloud Jaqaman, September 2009

%% Input

if nargin < 5
    disp('--plotMergeSplitPositions2D: Incorrect number of input arguments!');
    return
end

if nargin < 6 || isempty(figureName)
    figureName = [];
end

if nargin < 7 || isempty(image)
    image = [];
end

if nargin < 8 || isempty(mergesInfoSpaceExt)
    mergesInfoSpaceExt = repmat(mergesInfoSpace,[1 1 2]);
    mergesInfoSpaceExt(:) = NaN;
end

if nargin < 9 || isempty(splitsInfoSpaceExt)
    splitsInfoSpaceExt = repmat(splitsInfoSpace,[1 1 2]);
    splitsInfoSpaceExt(:) = NaN;
end

%get number of tracks
numTracks = length(tracks);

%% Pre-processing

%go over all tracks and calculate their center positions
%also determine the minimum and maximum coordinates
centerPosAllTracks = zeros(numTracks,2);
[minXCoord,minYCoord,maxXCoord,maxYCoord] = deal(0);
for iTrack = 1 : numTracks
    
    coordAmpCG = tracks(iTrack).tracksCoordAmpCG;
    xCoord = coordAmpCG(:,1:8:end);
    yCoord = coordAmpCG(:,2:8:end);
    
    centerPosAllTracks(iTrack,:) = [nanmean(xCoord(:)) nanmean(yCoord(:))];
    
    minXCoord = min(min(xCoord(:)),minXCoord);
    maxXCoord = max(max(xCoord(:)),maxXCoord);
    minYCoord = min(min(yCoord(:)),minYCoord);
    maxYCoord = max(max(yCoord(:)),maxYCoord);
    
end

%% Plotting

%make new figure
if isempty(figureName)
    figure
else
    figure('Name',figureName)
end

%plot the image to overlay the merge and split positions on, if given
if ~isempty(image)
    imshow(image,[]);
else
    imshow(ones(ceil(maxYCoord),ceil(maxXCoord)),[]);
end
hold on

%set figure axes limits
axis([minXCoord maxXCoord minYCoord maxYCoord]);

%show coordinates on axes
axH = gca;
set(axH,'visible','on');

%label axes
xlabel('x-coordinate (pixels)');
ylabel('y-coordinate (pixels)');

%first plot the center points of the tracks that do not undergo
%merges or splits

%find which tracks do not undergo merges or splits
indxTracksMS = unique([mergesInfo(:,1); splitsInfo(:,1)]);
indxTracksNoMS = setdiff(1:numTracks,indxTracksMS);

%plot the no M/S center positions
plot(centerPosAllTracks(indxTracksNoMS,1),centerPosAllTracks(indxTracksNoMS,2),'.','Color',[0.7 0.7 0.7])

%then plot the coordinates of splitting events

%get coordinates of splits
xCoord = splitsInfoSpace(:,1:2:end);
xCoord = xCoord(:);
yCoord = splitsInfoSpace(:,2:2:end);
yCoord = yCoord(:);
xCoord1 = splitsInfoSpaceExt(:,1:2:end,1);
xCoord1 = xCoord1(:);
yCoord1 = splitsInfoSpaceExt(:,2:2:end,1);
yCoord1 = yCoord1(:);
xCoord2 = splitsInfoSpaceExt(:,1:2:end,2);
xCoord2 = xCoord2(:);
yCoord2 = splitsInfoSpaceExt(:,2:2:end,2);
yCoord2 = yCoord2(:);

%keep only non-zero entries
indxKeep = find( xCoord~=0 & yCoord~=0 );
xCoord = xCoord(indxKeep);
yCoord = yCoord(indxKeep);
xCoord1 = xCoord1(indxKeep);
yCoord1 = yCoord1(indxKeep);
xCoord2 = xCoord2(indxKeep);
yCoord2 = yCoord2(indxKeep);

%format data for plotting
xCoord01 = [xCoord xCoord1]';
yCoord01 = [yCoord yCoord1]';
xCoord02 = [xCoord xCoord2]';
yCoord02 = [yCoord yCoord2]';

%plot the splits
lineWithGaps(xCoord01,yCoord01,[],'Color','k','Marker','.');
lineWithGaps(xCoord02,yCoord02,[],'Color','k','Marker','.','LineStyle',':');

%then plot the coordinates of merging events

%get coordinates of merges
xCoord = mergesInfoSpace(:,1:2:end);
xCoord = xCoord(:);
yCoord = mergesInfoSpace(:,2:2:end);
yCoord = yCoord(:);
xCoord1 = mergesInfoSpaceExt(:,1:2:end,1);
xCoord1 = xCoord1(:);
yCoord1 = mergesInfoSpaceExt(:,2:2:end,1);
yCoord1 = yCoord1(:);
xCoord2 = mergesInfoSpaceExt(:,1:2:end,2);
xCoord2 = xCoord2(:);
yCoord2 = mergesInfoSpaceExt(:,2:2:end,2);
yCoord2 = yCoord2(:);

%keep only non-zero entries
indxKeep = find( xCoord~=0 & yCoord~=0 );
xCoord = xCoord(indxKeep);
yCoord = yCoord(indxKeep);
xCoord1 = xCoord1(indxKeep);
yCoord1 = yCoord1(indxKeep);
xCoord2 = xCoord2(indxKeep);
yCoord2 = yCoord2(indxKeep);

%format data for plotting
xCoord01 = [xCoord xCoord1]';
yCoord01 = [yCoord yCoord1]';
xCoord02 = [xCoord xCoord2]';
yCoord02 = [yCoord yCoord2]';

%plot the merges
lineWithGaps(xCoord01,yCoord01,[],'Color','r','Marker','.');
lineWithGaps(xCoord02,yCoord02,[],'Color','r','Marker','.','LineStyle',':');

%make figure legend
legend('no merges/splits','splits1','splits2','merges1','merges2')


