function plotMergeSplitPositions2D(tracks,mergesInfo,splitsInfo,...
    mergesInfoSpace,splitsInfoSpace,figureName,image)
%PLOTMERGESPLITPOSITIONS2D generates a spatial map of merge and split locations
%
%SYNPOSIS plotMergeSplitPositions2D(tracks,mergesInfo,splitsInfo,...
%    mergesInfoSpace,splitsInfoSpace,figureName,image)
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
%                        storing the (x,y,[z])-coordinates of each merge.
%                        Every row corresponds to the same row in
%                        mergesInfo. Every merge gets 2 (in 2D) or 3 (in
%                        3D) columns for x, y and z (if 3D).
%       splitsInfoSpace: 2D array that is a continuation of splitsInfoTime,
%                        storing the (x,y,[z])-coordinates of each split.
%                        Every row corresponds to the same row in
%                        splitsInfo. Every split gets 2 (in 2D) or 3 (in
%                        3D) columns for x, y and z (if 3D).
%       figureName     : Figure name.
%                        Optional. Default: None.
%       image          : Image to overlay spatial map on.
%                        Optional. Default: [].
%
%REMARKS Merges are plotted in blue, splits in red, while trajectories
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

%get number of tracks
numTracks = length(tracks);

%% Pre-processing

% %convert tracks from structure to matrix format
% tracksMat = convStruct2MatIgnoreMS(tracks);
% 
% %find minimum and maximum coordinates for plotting
% xCoord = tracksMat(:,1:8:end);
% minXCoord = min(floor(min(xCoord(:))),0);
% maxXCoord =  ceil(max(xCoord(:)));
% yCoord = tracksMat(:,2:8:end);
% minYCoord = min(floor(min(yCoord(:))),0);
% maxYCoord =  ceil(max(yCoord(:)));

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
% numTracksNoMS = length(indxTracksNoMS);
% 
% %go over these tracks and get their center positions
% centerPosNoMS = zeros(numTracksNoMS,2);
% for iTrack = 1 : numTracksNoMS
%     coordAmpCG = tracks(indxTracksNoMS(iTrack)).tracksCoordAmpCG;
%     xCoord = coordAmpCG(:,1:8:end);
%     yCoord = coordAmpCG(:,2:8:end);
%     centerPosNoMS(iTrack,:) = [nanmean(xCoord(:)) nanmean(yCoord(:))];
% end
% 
% %plot the no M/S center positions
% plot(centerPosNoMS(:,1),centerPosNoMS(:,2),'.','Color',[0.7 0.7 0.7])

%plot the no M/S center positions
plot(centerPosAllTracks(indxTracksNoMS,1),centerPosAllTracks(indxTracksNoMS,2),'.','Color',[0.7 0.7 0.7])

%then plot the coordinates of splitting events

%get coordinates of splits
xCoord = splitsInfoSpace(:,1:2:end);
xCoord = xCoord(:);
yCoord = splitsInfoSpace(:,2:2:end);
yCoord = yCoord(:);

%keep only non-zero entries
indxKeep = find( xCoord~=0 & yCoord~=0 );
xCoord = xCoord(indxKeep);
yCoord = yCoord(indxKeep);

%plot the split locations
plot(xCoord,yCoord,'k+')

%then plot the coordinates of merging events

%get coordinates of merges
xCoord = mergesInfoSpace(:,1:2:end);
xCoord = xCoord(:);
yCoord = mergesInfoSpace(:,2:2:end);
yCoord = yCoord(:);

%keep only non-zero entries
indxKeep = find( xCoord~=0 & yCoord~=0 );
xCoord = xCoord(indxKeep);
yCoord = yCoord(indxKeep);

%plot the merge locations
plot(xCoord,yCoord,'rx')
%make figure legend
legend('no merges/splits','splits','merges')


