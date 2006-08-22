function plotTracks2D(trackedFeatureInfo,timeRange,colorTime,indicateSE,...
    newFigure,image)
%PLOTTRACKS2D plots a group of tracks in 2D and allows user to click on them and extract track information
%
%SYNOPSIS plotTracks2D(trackedFeatureInfo,timeRange,colorTime,indicateSE,...
%    newFigure,image)
%
%INPUT  trackedFeatureInfo: Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%       timeRange         : 2-element row vector indicating time range to plot. 
%                           Optional. Default: whole movie.
%       colorTime         : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                           blue, red, etc.
%                           Optional. Default: 'k'.
%       indicateSE        : 1 if track starts and ends are to be indicated
%                           with circles and squares, respectively; 0
%                           otherwise. Optional. Default: 1.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image.
%
%OUTPUT no output variables, just the plot
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--plotTracks2D: Incorrect number of input arguments!');
    return
end

%get number of tracks and number of time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/8;

errFlag = 0;

%check whether a time range for plotting was input
if nargin < 2 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracks2D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether colorTime was input
if nargin < 3 || isempty(colorTime)
    colorTime = 'k';
end

%check whether indicateSE was input
if nargin < 4 || isempty(indicateSE)
    indicateSE = 1;
else
    if indicateSE ~= 0 && indicateSE ~= 1
        disp('plotTracks2D: indicateSE should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 5 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1
        disp('--plotTracks2D: newFigure should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether user supplied an image
if nargin < 6 || isempty(image)
    image = [];
end

%exit if there are problem in input variables
if errFlag
    disp('--plotTracks2D: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the number of time points to be plotted
numTimePlot = timeRange(2) - timeRange(1) + 1;

%get the x,y-coordinates of features in all tracks
tracksX = trackedFeatureInfo(:,1:8:end)';
tracksY = trackedFeatureInfo(:,2:8:end)';

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

%if the user wants to plot in a new figure window
if newFigure
    
    %open new figure window
    figure

    if ~isempty(image) %if user supplied an image
        imshow(image,[]); %plot the image
    else %if user did not supply an image
        imshow(ones(ceil(max(tracksY(:))),ceil(max(tracksX(:)))),[]); %plot an empty image
    end
    
    %show coordinates on axes
    ah = gca;
    set(ah,'visible','on');

    %label axes
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');

end

%hold on figure
hold on

if colorTime == '1' %if user wants to color-code time

    %plot tracks ignoring missing points as a dotted black line
    for i=1:numTracks
        obsAvail = find(~isnan(tracksXP(:,i)));
        plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:')
    end

    %get the fraction of each color in each time interval to be plotted
    numTimePlotOver2 = ceil((numTimePlot-1)/2); %needed to change blue color over time
    redVariation = (0:numTimePlot-2)'/(numTimePlot-2);
    greenVariation = (numTimePlot-2:-1:0)'/(numTimePlot-2);
    blueVariation = [(0:numTimePlotOver2-1)'/(numTimePlotOver2-1);...
        (numTimePlot-numTimePlotOver2-2:-1:0)'/(numTimePlot-numTimePlotOver2-1)];

    %get the overall color per time interval
    colorOverTime = [redVariation greenVariation blueVariation];

    %overlay tracks with color coding wherever a feature has been detected
    for i=1:numTimePlot-1
        plot(tracksXP(i:i+1,:),tracksYP(i:i+1,:),'color',colorOverTime(i,:));
    end

else
    
    %plot tracks ignoring missing points with the line color indicated
    for i=1:numTracks
        obsAvail = find(~isnan(tracksXP(:,i)));
        plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),colorTime)
    end

end %(if colorTime == '1' ... else ...) 
    
if indicateSE %if user wants to indicate starts and ends

    %find the beginning and end of each track
    for i=numTracks:-1:1
        timePoint = find(~isnan(tracksX(:,i)));
        startInfo(i,:) = [tracksX(timePoint(1),i) ...
            tracksY(timePoint(1),i) timePoint(1)];
        endInfo(i,:) = [tracksX(timePoint(end),i) ...
            tracksY(timePoint(end),i) timePoint(end)];
    end

    %place circles at track starts and squares at track ends if they happen to
    %be in the plotting region of interest
    if colorTime == '1'
        indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
        plot(startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o');
        indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
        plot(endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square');
    else
        indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
        plot(startInfo(indx,1),startInfo(indx,2),colorTime,...
            'LineStyle','none','marker','o');
        indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
        plot(endInfo(indx,1),endInfo(indx,2),colorTime,...
            'LineStyle','none','marker','square');
    end

end

%ask the user whether to click on figure and get frame information
userEntry = input('select points in figure? y/n ','s');

while strcmp(userEntry,'y')

    %let the user choose the points of interest
    [x,y] = getpts;

    %find the time points of the indicated points
    for i=1:length(x)
        distTrack2Point = (tracksXP-x(i)).^2+(tracksYP-y(i)).^2;
        [frameChosen,trackChosen] = find(distTrack2Point==min(distTrack2Point(:)));
        for j=1:length(trackChosen)
            disp(['Track: ' num2str(trackChosen(j)) ...
                '   Frame: ' num2str(frameChosen(j)+timeRange(1)-1) ...
                '   Coordinates: ' num2str(tracksXP(frameChosen(j),trackChosen(j))) ...
                ' ' num2str(tracksYP(frameChosen(j),trackChosen(j)))  ]);
        end
    end
        
    %ask the user again whether to click on figure and get frame information
    userEntry = input('select points again? y/n ','s');

end

%%%%% ~~ the end ~~ %%%%%

