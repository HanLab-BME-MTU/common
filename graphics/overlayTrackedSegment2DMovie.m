function overlayTrackedSegment2DMovie(tracksFinal,segmentParams,startend,...
    dragtailLength,inputMovieInfo,saveMovieInfo,classifyGaps,highlightES,...
    imageRange,onlyTracks,classifyLft,diffAnalysisRes,colorTracks)
%OVERLAYTRACKEDSEGMENT2DMOVIE Overlays tracked segment 2D obtained via
%trackCloseGapsKalman on movies with variable color-coding
%
%SYNPOSIS overlayTrackedSegment2DMovie(tracksFinal,segmentParams,startend,...
%    dragtailLength,inputMovieInfo, saveMovieInfo,classifyGaps,highlightES,...
%    imageRange,onlyTracks,classifyLft,diffAnalysisRes,colorTracks)
%
%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       segmentsInfo  : Cell array of size equal to the number of frames.
%                       Each cell contains the parameters of detection
%                       2D segments, i.e. an array of size Nx5, where N is
%                       the number of detected segments at a given frame
%                       and the number of columns corresponds to the xy
%                       position, amplitude, length and orientation of a
%                       segment.
%       startend      : Row vector indicating first and last frame to
%                       include in movie. Format: [startframe endframe].
%                       Optional. Default: [(first frame with tracks) (last frame with tracks)]
%       dragtailLength: Length of drag tail (in frames).
%                       Optional. Default: 10 frames.
%                       ** If dragtailLength = 0, then no dragtail.
%                       ** To show tracks from their beginning to their end,
%                       set dragtailLength to any value longer than the
%                       movie.
%                       ** To show tracks statically while features dance
%                       on them, use -1.
%                       ** To show tracks from their beginning to their
%                       end, and to retain tracks even after the particle
%                       disappears, use -2.
%       inputMovieInfo: structure providing the location of raw images to
%                       be overlaid.
%           .dir          : Directory where the raw images are located.
%           .filename     : filename of the very first image in the movie.
%       saveMovieInfo : empty struct if no saving is requested.
%           .dir          : Directory where the movie should be saved.
%           .filename     : Name of file where results should be saved.
%       classifyGaps  : 1 to classify gaps as "good" and "bad", depending
%                       on their length relative to the legnths of the
%                       segments they connect, 0 otherwise.
%                       Optional. Default: 1.
%       highlightES   : 1 to highlight track ends and starts, 0 otherwise.
%                       Optional. Default: 1.
%       imageRange    : Image region to make movie out of, in the form:
%                       [min pixel X, max pixel X; min pixel Y, max pixel Y].
%                       Optional. Default: Whole image.
%       onlyTracks    : 1 to show only tracks without any symbols showing
%                       detections, closed gaps, merges and splits; 0 to
%                       show symbols on top of tracks.
%                       Optional. Default: 0.
%       classifyLft   : 1 to classify objects based on (1) whether they
%                       exist throughout the whole movie, (2) whether they
%                       appear OR disappear, and (3) whether they appear
%                       AND disappear; 0 otherwise.
%                       Optional. Default: 1.
%       diffAnalysisRes: Diffusion analysis results (either output of
%                       trackDiffusionAnalysis1 or trackTransientDiffusionAnalysis2).
%                       Needed if tracks/track segments are to be
%                       colored based on their diffusion classification.
%                       With this option, classifyGaps, highlightES and
%                       classifyLft are force-set to zero, regardless of input.
%                       Optional. Default: None.
%       colorTracks   : 1 to color tracks by rotating through 7 different
%                       colors, 0 otherwise. With this option,
%                       classifyGaps, highlightES and classifyLft are
%                       force-set to zero, regardless of input.
%                       Option ignored if diffAnalysisRes is supplied.
%                       Optional. Default: 0.
%
%REMARKS Color-coding:
%        ** Without diffusion classification, all tracks have a neutral
%        color, while objects are color coded in the following way:
%               * Detected object just after appearance: Green circle.
%               * Detected object just before disappearance: Yellow
%                 circle.
%               * Detected object in middle of trajectory that spans
%                 whole movie: White circle.
%               * Detected object in middle of trajectory that appears OR
%                 disappears within movie: Magenta circle.
%               * Detected object in middle of trajectory that appears AND
%                 disappears within movie: Red circle.
%               * Gap that is short than both segments it connects: Cyan
%                 star.
%               * Gap that is longer than at least one ofthe segments it
%                 connects: Blue star.
%               * Object before and after splitting: Green diamond.
%               * OBject before and after merging: Yellow diamond.
%           When classifyGaps = 0, all gaps are cyan.
%           When highlighES = 0, no green and yellow circles.
%           When classifyLft = 0, all objets in middle of trajectory are white.
%
%       ** With diffusion classification, all objects and gaps have neutral
%       color (merges and splits are diamonds), while tracks and track
%       segments are color-coded in the following way:
%               * Type 1: Linear + 1D confined diffusion: Orange.
%               * Type 2: Linear + 1D normal diffusion: Red.
%               * Type 3: Linear + 1D super diffusion: Green.
%               * Type 4: Linear + too short for 1D classification: Yellow.
%               * Type 5: Random/Unclassified + 2D confined diffusion: Blue.
%               * Type 6: Random/Unclassified + 2D normal diffusion: Cyan.
%               * Type 7: Random/Unclassified + 2D super diffusion: Magenta.
%               * Type 8: Random + too short for 2D classification: Purple.
%               * Type 0: Too short for any analysis: Light pink.
%
%NOTE this file has been adapted from overlayTracksMovieNew.m. Some changes
%have been made:
% - not every parameter of the tracked segments are stored in tracksFinal.
%   We need therefore to provide segmentParams which contains the full set
%   of segment parameters. To access parameters, we use tracksFeatIndxCG
%   instead of tracksCoordAmpCG.
% - parameters filterSigma, showRaw and intensityScale have been removed.
% - parameter saveMovie and movieName have been bundled into a structure
% saveMovieInfo containing the directory and the filename of the saved
% movie. if the structure is empty, the movie is not saved. There is no
% more a default filename and default location in case the movie wants to
% be saved. This new parameter avoids to save the movie into directory
% where raw images are located.
% - inputMovieInfo is a structure containing the directory and the first
% name of the raw images to be overlaid. In case inputMovieInfo is not
% provided, the user will be asked to specify first image of the movie,
% complying with the initial behavious or overlayFeaturesMovie.m. This new
% parameter allows to batch movie creation.
%
%Khuloud Jaqaman, August 2007
%Adapted by Sylvain Berlemont, 2010

%% input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--overlayTrackedSegment2DMovie: Incorrect number of input arguments!');
    return
end

%get first and last frames where there are tracks
allEvents = vertcat(tracksFinal.seqOfEvents);
tracksFirstFrame = min(allEvents(:,1));
tracksLastFrame = max(allEvents(:,1));

%check whether segment parameters are provided to the last frame
if numel(segmentParams) < tracksLastFrame
    disp('--overSegment2DMovie: Not enough segment parameters.');
    return;
end

%check whether a inputMovieInfo has been provided
if nargin < 5 || isempty(inputMovieInfo)
    %ask user for images
    [fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image, even if not to be plotted');

else
    if ~isfield(inputMovieInfo, 'dir') || ~exist(inputMovieInfo.dir, 'dir')
        disp('--overSegment2DMovie: Invalid input directory.');
        return;
    end
    
    if ~isfield(inputMovieInfo, 'filename') || ~ischar(inputMovieInfo.filename)
        disp('--overSegment2DMovie: Invalid input filename.');
        return;
    end
    
    fName = inputMovieInfo.filename;
    dirName = inputMovieInfo.dir;
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames(fullfile(dirName,fName));
    numFiles = length(outFileList);
    
    %determine which frames the files correspond to, and generate the inverse map
    %indicate missing frames with a zero
    frame2fileMap = zeros(numFiles,1);
    for iFile = 1 : numFiles
        [~,~,frameNumStr] = getFilenameBody(outFileList{iFile});
        frameNum = str2double(frameNumStr);
        frame2fileMap(frameNum) = iFile;
    end
    
    %assign as number of frames the last frame number observed
    numFrames = frameNum;
    
    %read first image to get image size
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    
else %else, exit 
    disp('--overlayTrackedSegment2DMovie: Bad file selection');
    return
    
end

%check startend and assign default if necessary
if nargin < 3 || isempty(startend)
    startend = [tracksFirstFrame tracksLastFrame];
else
    tracksFirstFrame = min(tracksFirstFrame,startend(1));
    tracksLastFrame = max(tracksLastFrame,startend(2));
end

%keep only the frames of interest
outFileList = outFileList(frame2fileMap(startend(1)):frame2fileMap(startend(2)));
frame2fileMap = frame2fileMap(startend(1):startend(2));
indxNotZero = find(frame2fileMap~=0);
frame2fileMap(indxNotZero) = frame2fileMap(indxNotZero) - frame2fileMap(indxNotZero(1)) + 1;

%check dragtailLength and assign default if not necessary
if nargin < 4 || isempty(dragtailLength)
    dragtailLength = 10;
end

%check whether to save movie
if nargin < 6 || isempty(saveMovieInfo)
    saveMovie = false;
else
    if ~isfield(saveMovieInfo, 'dir') || ~exist(saveMovieInfo.dir, 'dir')
        disp('--overSegment2DMovie: Invalid output directory.');
        return;
    end
    
    if ~isfield(saveMovieInfo, 'filename') || ~ischar(saveMovieInfo.filename)
        disp('--overSegment2DMovie: Invalid output filename.');
        return;
    end
    
    [~,~,~,ext] = getFilenameBody(saveMovieInfo.filename);
    
    if ~strcmp(ext, '.mov')
        disp('--overSegment2DMovie: Invalid output filename format (.mov is required).');
        return;
    end
    
    saveMovie = true;
end

%check whether to color-code gaps
if nargin < 7 || isempty(classifyGaps)
    classifyGaps = 1;
end

%check whether to highligh track starts and ends
if nargin < 8 || isempty(highlightES)
    highlightES = 1;
end

%check whether an area of interest was input
if nargin < 9 || isempty(imageRange)
    imageRange = [1 isx; 1 isy];
end

%check whether to plot tracks only or also symbols
if nargin < 10 || isempty(onlyTracks)
    onlyTracks = 0;
end

%check whether to classify lifetime
if nargin < 11 || isempty(classifyLft)
    classifyLft = 1;
end

%check whether to color-code tracks based on diffusion classification
%check whether diffusion classification is for overall tracks or transient
if nargin < 12 || isempty(diffAnalysisRes)
    diffAnalysisRes = [];
    transDiffClass = 0;
else
    classifyGaps = 0;
    highlightES = 0;
    classifyLft = 0;
    if isfield(diffAnalysisRes,'segmentClass')
        transDiffClass = 1;
    else
        transDiffClass = 0;
    end
end

%check whether to color individual tracks
if nargin < 13 || isempty(colorTracks)
    colorTracks = 0;
else
    if ~isempty(diffAnalysisRes)
        colorTracks = 0;
    end
    if colorTracks == 1
        classifyGaps = 0;
        highlightES = 0;
        classifyLft = 0;
    end
end

%initialize QT movie if it is to be saved
if saveMovie
    evalString = ['MakeQTMovie start ''' fullfile(saveMovieInfo.dir,saveMovieInfo.filename) ''''];
    eval(evalString);
end

%% store track positions, get track status and point status

%get number of tracks
numTracks = length(tracksFinal);

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%give tracks status based on the frames they span:
%2: track exists throughout movie
%1: track exists either in first frame or in last frame
%0: track does not exist in both first frame and last frame
trackStatus  = (trackSEL(:,1) == tracksFirstFrame) + (trackSEL(:,2) == tracksLastFrame);

%give all tracks same classification if lifetime classification not
%requested
if classifyLft == 0
    trackStatus(:) = 2;
end

%get number of segments making each track
numSegments = zeros(numTracks,1);
for i = 1 : numTracks
    numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
end

%locate the row of the first segment of each compound track in the
%big matrices of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%find total number of segments in all tracks (i.e. number of rows in big
%matrices)
numSegmentsTracks = trackStartRow(end)+numSegments(end)-1;

%construct a matrix indicating point status in big matrices:
%-2: bad gap (gap length > either segment length on its sides)
%-1: good gap (gap length < both segment lengths on its sides)
%0 : before track start or after track end
%1 : detected feature in the middle of a track of trackStatus = 0
%2 : detected feature in the middle of a track of trackStatus = 1
%3 : detected feature in the middle of a track of trackStatus = 2
%4 : detected feature just after a birth
%5 : detected feature just before a death
%6 : detected feature just after a split
%7 : detected feature just before a merge
pointStatus = zeros(numSegmentsTracks,numFrames);

%put all tracks together in one big matrix
%put the x-coordinates in one matrix and the y-coordinates in another
%indicate the status of each point

xCoordMatAll = NaN*ones(numSegmentsTracks,numFrames);
yCoordMatAll = xCoordMatAll;
lengthMatAll = xCoordMatAll;
angleMatAll = xCoordMatAll;

for iTrack = 1 : numTracks
    
    %get track start and end times
    startTime = trackSEL(iTrack,1);
    endTime   = trackSEL(iTrack,2);
    
    %get feature indices of current track
    tracksFeatIndxCG = tracksFinal(iTrack).tracksFeatIndxCG;
    
    %store x-coordinates
    xCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = arrayfun(@(iFrame) ...
        segmentParams{iFrame}(tracksFeatIndxCG(:,iFrame-startTime+1), 1),...
        startTime:endTime);
    
    %store y-coordinates
    yCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = arrayfun(@(iFrame) ...
        segmentParams{iFrame}(tracksFeatIndxCG(:,iFrame-startTime+1), 2),...
        startTime:endTime);

    %store lengths
    lengthMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = arrayfun(@(iFrame) ...
        segmentParams{iFrame}(tracksFeatIndxCG(:,iFrame-startTime+1), 4),...
        startTime:endTime);
    
    %store orientations
    angleMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = arrayfun(@(iFrame) ...
        segmentParams{iFrame}(tracksFeatIndxCG(:,iFrame-startTime+1), 5),...
        startTime:endTime);
    
    %assign point status for features in the middle of the track
    pointStatus(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = trackStatus(iTrack) + 1;
    
    %get sequence of events of track
    seqOfEvents = tracksFinal(iTrack).seqOfEvents;
    
    if highlightES
        
        %assign point status for features just after a birth
        points2consider = find(seqOfEvents(:,2)==1 & isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1)~=tracksFirstFrame)';
        for iPoint = points2consider
            pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                seqOfEvents(iPoint,1)) = 4;
        end
        
        %assign point status for features just before a death
        points2consider = find(seqOfEvents(:,2)==2 & isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1)~=tracksLastFrame)';
        for iPoint = points2consider
            pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                seqOfEvents(iPoint,1)) = 5;
        end
        
    end
    
    %assign point status for features just after and before a split
    %also, in the frame just before splitting, give the
    %splitting track the position of the track it split from
    points2consider = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
        lengthMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = lengthMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
        angleMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = angleMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
    end
    
    %assign point status for features just before and after a merge
    %also, in the frame just after merging, give the
    %merging track the position of the track it merged from
    points2consider = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
        lengthMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = lengthMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
        angleMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = angleMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
    end
end %(for iTrack = 1 : numTracks)

%find gaps in tracks
gapInfo = findTrackGaps(tracksFinal);

%for gaps, assign the position as that of the feature before the gap
%also, assign a point status of -1
for iGap = 1 : size(gapInfo,1)
    
    iTrack = gapInfo(iGap,1);
    iSegment = gapInfo(iGap,2);
    iStart = gapInfo(iGap,3);
    gapLength = gapInfo(iGap,4);
    if classifyGaps
        if gapInfo(iGap,5) <= 1 && gapInfo(iGap,6) <= 1
            gapType = -1;
        else
            gapType = -2;
        end
    else
        gapType = -1;
    end
    
    xCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        xCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %x-coordinates
    
    yCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        yCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %y-coordinates
    
    lengthMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        lengthMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %lengths
    
    angleMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        angleMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %orientations
    
    pointStatus(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = gapType; %point status
end

%retain in the big matrices only the frames of interest
xCoordMatAll = xCoordMatAll(:,startend(1):startend(2));
yCoordMatAll = yCoordMatAll(:,startend(1):startend(2));
lengthMatAll = lengthMatAll(:,startend(1):startend(2));
angleMatAll = angleMatAll(:,startend(1):startend(2));

pointStatus = pointStatus(:,startend(1):startend(2));

%% divide tracks based on diffusion analysis or just into groups to be colored separately

%if tracks are to be individually colored
if colorTracks
    
    %divide tracks among 9 matrices that will get their own colors
    
    %x-coordinates ...
    [xCoordMatAll0,xCoordMatAll1,xCoordMatAll2,xCoordMatAll3,xCoordMatAll4,...
        xCoordMatAll5,xCoordMatAll6,xCoordMatAll7,xCoordMatAll8] = deal(NaN(size(xCoordMatAll)));
    xCoordMatAll0(1:9:end,:) = xCoordMatAll(1:9:end,:);
    xCoordMatAll1(2:9:end,:) = xCoordMatAll(2:9:end,:);
    xCoordMatAll2(3:9:end,:) = xCoordMatAll(3:9:end,:);
    xCoordMatAll3(4:9:end,:) = xCoordMatAll(4:9:end,:);
    xCoordMatAll4(5:9:end,:) = xCoordMatAll(5:9:end,:);
    xCoordMatAll5(6:9:end,:) = xCoordMatAll(6:9:end,:);
    xCoordMatAll6(7:9:end,:) = xCoordMatAll(7:9:end,:);
    xCoordMatAll7(8:9:end,:) = xCoordMatAll(8:9:end,:);
    xCoordMatAll8(9:9:end,:) = xCoordMatAll(9:9:end,:);
    
    %y-coordinates ...
    [yCoordMatAll0,yCoordMatAll1,yCoordMatAll2,yCoordMatAll3,yCoordMatAll4,...
        yCoordMatAll5,yCoordMatAll6,yCoordMatAll7,yCoordMatAll8] = deal(NaN(size(yCoordMatAll)));
    yCoordMatAll0(1:9:end,:) = yCoordMatAll(1:9:end,:);
    yCoordMatAll1(2:9:end,:) = yCoordMatAll(2:9:end,:);
    yCoordMatAll2(3:9:end,:) = yCoordMatAll(3:9:end,:);
    yCoordMatAll3(4:9:end,:) = yCoordMatAll(4:9:end,:);
    yCoordMatAll4(5:9:end,:) = yCoordMatAll(5:9:end,:);
    yCoordMatAll5(6:9:end,:) = yCoordMatAll(6:9:end,:);
    yCoordMatAll6(7:9:end,:) = yCoordMatAll(7:9:end,:);
    yCoordMatAll7(8:9:end,:) = yCoordMatAll(8:9:end,:);
    yCoordMatAll8(9:9:end,:) = yCoordMatAll(9:9:end,:);
else %otherwise
    
    %copy all coordinates into new variables
    %only the "0" variables will be used when there is no diffusion
    %classification
    [xCoordMatAll0,xCoordMatAll1,xCoordMatAll2,xCoordMatAll3,xCoordMatAll4,...
        xCoordMatAll5,xCoordMatAll6,xCoordMatAll7,xCoordMatAll8] = deal(xCoordMatAll);
    [yCoordMatAll0,yCoordMatAll1,yCoordMatAll2,yCoordMatAll3,yCoordMatAll4,...
        yCoordMatAll5,yCoordMatAll6,yCoordMatAll7,yCoordMatAll8] = deal(yCoordMatAll);
end

%if there are diffusion analysis results ...
if ~isempty(diffAnalysisRes)
    
    if transDiffClass %if transient diffusion classification ...
        
        %put all track segment classifications into one array (note that
        %here we're talking about merging and splitting segments)
        trackSegmentClass = vertcat(diffAnalysisRes.segmentClass);
        
        %get number of time points in movie
        numTimePoints = size(xCoordMatAll,2);
        
        %go over all track segments and classify their points
        for iTrackSegment = 1 : numSegmentsTracks
            
            %get current track's classification
            trackClassCurrent = trackSegmentClass(iTrackSegment).momentScalingSpectrum(:,1:3);
            
            %map the transient classification numbers into the same numbers
            %as the whole track classification
            trackClassCol3 = trackClassCurrent(:,3);
            trackClassCurrent(trackClassCol3==1,3) = 5; %random/unclassified + 2D confined (blue)
            trackClassCurrent(trackClassCol3==2,3) = 6; %random/unclassified + 2D normal (cyan)
            trackClassCurrent(trackClassCol3==3,3) = 7; %random/unclassified + 2D super (magenta)
            trackClassCurrent(isnan(trackClassCol3),3) = 0; %completely unclassified
            
            %extract and store segments classified as 5
            classCurrent = find(trackClassCurrent(:,3)==5); %find points classified as 5
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 5
            xCoordMatAll5(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 5
            yCoordMatAll5(iTrackSegment,pointsNotCurrent) = NaN;
            
            %extract and store segments classified as 6
            classCurrent = find(trackClassCurrent(:,3)==6); %find points classified as 6
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 6
            xCoordMatAll6(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 6
            yCoordMatAll6(iTrackSegment,pointsNotCurrent) = NaN;
            
            %extract and store segments classified as 7
            classCurrent = find(trackClassCurrent(:,3)==7); %find points classified as 7
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 7
            xCoordMatAll7(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 7
            yCoordMatAll7(iTrackSegment,pointsNotCurrent) = NaN;
            
            %extract and store unclassified segments
            classCurrent = find(trackClassCurrent(:,3)==0); %find points classified as 0
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 0
            xCoordMatAll0(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 0
            yCoordMatAll0(iTrackSegment,pointsNotCurrent) = NaN;       
        end
        
    else %if whole track diffusion classification ...
        
        %get track segment types from diffusion analysis
        trackSegmentType = vertcat(diffAnalysisRes.classification);
        
        %assign classes
        trackClass = zeros(numSegmentsTracks,1); %initialize with the indicator for undetermined
        trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 1) = 1; %linear + 1D confined (orange)
        trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 2) = 2; %linear + 1D normal (red)
        trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 3) = 3; %linear + 1D super (green)
        trackClass(trackSegmentType(:,1) == 1 & isnan(trackSegmentType(:,3))) = 4; %linear + too short (yellow)
        trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 1) = 5; %random/unclassified + 2D confined (blue)
        trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 2) = 6; %random/unclassified + 2D normal (cyan)
        trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 3) = 7; %random/unclassified + 2D super (magenta)
        trackClass(trackSegmentType(:,1) == 0 & isnan(trackSegmentType(:,2))) = 8; %random + too short (purple)
        
        %extract the tracks/track segments of different classifications
        %x-coordinates
        xCoordMatAll0(trackClass~=0,:) = NaN;
        xCoordMatAll1(trackClass~=1,:) = NaN;
        xCoordMatAll2(trackClass~=2,:) = NaN;
        xCoordMatAll3(trackClass~=3,:) = NaN;
        xCoordMatAll4(trackClass~=4,:) = NaN;
        xCoordMatAll5(trackClass~=5,:) = NaN;
        xCoordMatAll6(trackClass~=6,:) = NaN;
        xCoordMatAll7(trackClass~=7,:) = NaN;
        xCoordMatAll8(trackClass~=8,:) = NaN;
        %y-coordinates
        yCoordMatAll0(trackClass~=0,:) = NaN;
        yCoordMatAll1(trackClass~=1,:) = NaN;
        yCoordMatAll2(trackClass~=2,:) = NaN;
        yCoordMatAll3(trackClass~=3,:) = NaN;
        yCoordMatAll4(trackClass~=4,:) = NaN;
        yCoordMatAll5(trackClass~=5,:) = NaN;
        yCoordMatAll6(trackClass~=6,:) = NaN;
        yCoordMatAll7(trackClass~=7,:) = NaN;
        yCoordMatAll8(trackClass~=8,:) = NaN;
        
    end %(if transDiffClass ... else ...)
    
end

%% make movie

%go over all specified frames
for iFrame = 1 : size(xCoordMatAll,2)
    
    if frame2fileMap(iFrame) ~= 0 %if frame exists 
        %read specified image
        imageStack = imread(outFileList{frame2fileMap(iFrame)});   
    else %otherwise
        %make empty frame
        imageStack = zeros(isx,isy);    
    end
    
    %plot image in current frame and show frame number
    clf;
    
    imagesc(imageStack),colormap gray,axis image,axis off;
    hold on;
    textDeltaCoord = min(diff(imageRange,[],2))/20;
    text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
        textDeltaCoord,num2str(iFrame+startend(1)-1),'Color','white');
    
    %get tracks to plot
    plotOrNot = 0;
    if dragtailLength >= 0 %to plot tracks dynamically
        
        if iFrame > 1 || onlyTracks
            dragTailStart = max(iFrame-dragtailLength,1);
            indx2keep = find(pointStatus(:,iFrame)~=0);
            xCoord2plot0 = (xCoordMatAll0(indx2keep,dragTailStart:iFrame))';
            yCoord2plot0 = (yCoordMatAll0(indx2keep,dragTailStart:iFrame))';
            xCoord2plot1 = (xCoordMatAll1(indx2keep,dragTailStart:iFrame))';
            yCoord2plot1 = (yCoordMatAll1(indx2keep,dragTailStart:iFrame))';
            xCoord2plot2 = (xCoordMatAll2(indx2keep,dragTailStart:iFrame))';
            yCoord2plot2 = (yCoordMatAll2(indx2keep,dragTailStart:iFrame))';
            xCoord2plot3 = (xCoordMatAll3(indx2keep,dragTailStart:iFrame))';
            yCoord2plot3 = (yCoordMatAll3(indx2keep,dragTailStart:iFrame))';
            xCoord2plot4 = (xCoordMatAll4(indx2keep,dragTailStart:iFrame))';
            yCoord2plot4 = (yCoordMatAll4(indx2keep,dragTailStart:iFrame))';
            xCoord2plot5 = (xCoordMatAll5(indx2keep,dragTailStart:iFrame))';
            yCoord2plot5 = (yCoordMatAll5(indx2keep,dragTailStart:iFrame))';
            xCoord2plot6 = (xCoordMatAll6(indx2keep,dragTailStart:iFrame))';
            yCoord2plot6 = (yCoordMatAll6(indx2keep,dragTailStart:iFrame))';
            xCoord2plot7 = (xCoordMatAll7(indx2keep,dragTailStart:iFrame))';
            yCoord2plot7 = (yCoordMatAll7(indx2keep,dragTailStart:iFrame))';
            xCoord2plot8 = (xCoordMatAll8(indx2keep,dragTailStart:iFrame))';
            yCoord2plot8 = (yCoordMatAll8(indx2keep,dragTailStart:iFrame))';
            plotOrNot = 1;
        end
        
    elseif dragtailLength == -1 %to plot tracks statically
        
        xCoord2plot0 = (xCoordMatAll0)';
        yCoord2plot0 = (yCoordMatAll0)';
        xCoord2plot1 = (xCoordMatAll1)';
        yCoord2plot1 = (yCoordMatAll1)';
        xCoord2plot2 = (xCoordMatAll2)';
        yCoord2plot2 = (yCoordMatAll2)';
        xCoord2plot3 = (xCoordMatAll3)';
        yCoord2plot3 = (yCoordMatAll3)';
        xCoord2plot4 = (xCoordMatAll4)';
        yCoord2plot4 = (yCoordMatAll4)';
        xCoord2plot5 = (xCoordMatAll5)';
        yCoord2plot5 = (yCoordMatAll5)';
        xCoord2plot6 = (xCoordMatAll6)';
        yCoord2plot6 = (yCoordMatAll6)';
        xCoord2plot7 = (xCoordMatAll7)';
        yCoord2plot7 = (yCoordMatAll7)';
        xCoord2plot8 = (xCoordMatAll8)';
        yCoord2plot8 = (yCoordMatAll8)';
        plotOrNot = 1;
        
    elseif dragtailLength == -2 %to plot tracks dynamically but keep them after they disappear
        
        if iFrame > 1 || onlyTracks
            xCoord2plot0 = xCoordMatAll0(:,1:iFrame)';
            yCoord2plot0 = yCoordMatAll0(:,1:iFrame)';
            xCoord2plot1 = xCoordMatAll1(:,1:iFrame)';
            yCoord2plot1 = yCoordMatAll1(:,1:iFrame)';
            xCoord2plot2 = xCoordMatAll2(:,1:iFrame)';
            yCoord2plot2 = yCoordMatAll2(:,1:iFrame)';
            xCoord2plot3 = xCoordMatAll3(:,1:iFrame)';
            yCoord2plot3 = yCoordMatAll3(:,1:iFrame)';
            xCoord2plot4 = xCoordMatAll4(:,1:iFrame)';
            yCoord2plot4 = yCoordMatAll4(:,1:iFrame)';
            xCoord2plot5 = xCoordMatAll5(:,1:iFrame)';
            yCoord2plot5 = yCoordMatAll5(:,1:iFrame)';
            xCoord2plot6 = xCoordMatAll6(:,1:iFrame)';
            yCoord2plot6 = yCoordMatAll6(:,1:iFrame)';
            xCoord2plot7 = xCoordMatAll7(:,1:iFrame)';
            yCoord2plot7 = yCoordMatAll7(:,1:iFrame)';
            xCoord2plot8 = xCoordMatAll8(:,1:iFrame)';
            yCoord2plot8 = yCoordMatAll8(:,1:iFrame)';
            plotOrNot = 1;
        end
        
    end
    
    %plot tracks
    if plotOrNot
        
        %plot basic tracks
        line([xCoord2plot0(1,:); xCoord2plot0],[yCoord2plot0(1,:); yCoord2plot0],...
            'Color',[1 0.7 0.7],'LineWidth',1); %light pink; the artificial repetition of the first line is for avoiding a mess in the first frame when tracks are not color-coded individually
        
        %color individual tracks randomly if requested
        if colorTracks == 1
            line(xCoord2plot1,yCoord2plot1,'Color',[1 0.7 0],'LineWidth',1); %orange
            line(xCoord2plot2,yCoord2plot2,'Color','r','LineWidth',1); %[1 0 0]
            line(xCoord2plot3,yCoord2plot3,'Color','g','LineWidth',1); %[0 1 0]
            line(xCoord2plot4,yCoord2plot4,'Color','y','LineWidth',1); %[1 1 0]
            line(xCoord2plot5,yCoord2plot5,'Color','b','LineWidth',1); %[0 0 1]
            line(xCoord2plot6,yCoord2plot6,'Color','c','LineWidth',1); %[0 1 1]
            line(xCoord2plot7,yCoord2plot7,'Color','m','LineWidth',1); %[1 0 1]
            line(xCoord2plot8,yCoord2plot8,'Color',[0.6 0 1],'LineWidth',1); %purple
        end
        
        %color-code dragtail based on diffusion analysis if supplied
        if ~isempty(diffAnalysisRes)
            line(xCoord2plot5,yCoord2plot5,'Color','b','LineWidth',1); %[0 0 1]
            line(xCoord2plot6,yCoord2plot6,'Color','c','LineWidth',1); %[0 1 1]
            line(xCoord2plot7,yCoord2plot7,'Color','m','LineWidth',1); %[1 0 1]
            if ~transDiffClass
                line(xCoord2plot1,yCoord2plot1,'Color',[1 0.7 0],'LineWidth',1); %orange
                line(xCoord2plot2,yCoord2plot2,'Color','r','LineWidth',1); %[1 0 0]
                line(xCoord2plot3,yCoord2plot3,'Color','g','LineWidth',1); %[0 1 0]
                line(xCoord2plot4,yCoord2plot4,'Color','y','LineWidth',1); %[1 1 0]
                line(xCoord2plot8,yCoord2plot8,'Color',[0.6 0 1],'LineWidth',1); %purple
            end
        end
        
    end
    
    %plot points (segment + gaps + merges + splits)
    if ~onlyTracks
        
        %blue stars: bad gaps
        points2plot = find(pointStatus(:,iFrame)==-2);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'b');
        line(xC,yC,'LineStyle','none','Color','b','Marker','*','MarkerSize',6);
        
        %cyan stars: good gaps
        points2plot = find(pointStatus(:,iFrame)==-1);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        if isempty(diffAnalysisRes)
            line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'c');
            line(xC,yC,'LineStyle','none','Color','c','Marker','*','MarkerSize',6);
        else
            line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'w');
            line(xC,yC,'LineStyle','none','Color','w','Marker','*','MarkerSize',6);
        end
        
        %red circles: detected feature in the middle of track with status 0
        points2plot = find(pointStatus(:,iFrame)==1);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);

        line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'r');
        line(xC,yC,'LineStyle','none','Color','r','Marker','o','MarkerSize',5);
        
        %magenta circles: detected feature in the middle of track with status 1
        points2plot = find(pointStatus(:,iFrame)==2);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'm');
        line(xC,yC,'LineStyle','none','Color','m','Marker','o','MarkerSize',5);
        
        %white circles: detected feature in the middle of track with status 2
        points2plot = find(pointStatus(:,iFrame)==3);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'w');
        line(xC,yC,'LineStyle','none','Color','w','Marker','o','MarkerSize',5);
        
        %green circles: detected feature just after birth
        points2plot = find(pointStatus(:,iFrame)==4);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'g');
        line(xC,yC,'LineStyle','none','Color','g','Marker','o','MarkerSize',5);
        
        %yellow circles: detected feature just before death
        points2plot = find(pointStatus(:,iFrame)==5);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'y');
        line(xC,yC,'LineStyle','none','Color','y','Marker','o','MarkerSize',5);
        
        %green diamonds: detected feature just before/after a split
        points2plot = find(pointStatus(:,iFrame)==6);
        xC = xCoordMatAll(points2plot,iFrame);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        if isempty(diffAnalysisRes)
            line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'g');
            line(xC,yC,'LineStyle','none','Color','g','Marker','d','MarkerSize',6);
        else
            line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'w');
            line(xC,yC,'LineStyle','none','Color','w','Marker','d','MarkerSize',6);
        end
        
        %yellow diamonds: detected feature just before/after a merge
        points2plot = find(pointStatus(:,iFrame)==7);
        yC = yCoordMatAll(points2plot,iFrame);
        l = lengthMatAll(points2plot,iFrame);
        t = angleMatAll(points2plot,iFrame);
        ct = cos(t);
        st = sin(t);
        
        if isempty(diffAnalysisRes)
            line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'y');
            line(xC,yC,'LineStyle','none','Color','y','Marker','d','MarkerSize',6);
        else
            line([xC - (l / 2) .* ct, xC + (l / 2) .* ct]', [yC - (l / 2) .* st, yC + (l / 2) .* st]', 'Color', 'w');
            line(xC,yC,'LineStyle','none','Color','w','Marker','d','MarkerSize',6);
        end
        
    end
    
    %add frame to movie if movie is saved
    if saveMovie
        MakeQTMovie addaxes
    end
    
end

%finish movie
if saveMovie==1
    MakeQTMovie finish
end
