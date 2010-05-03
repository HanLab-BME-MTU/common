function overlaySegment2DMovie(segmentParams,startend,saveMovieInfo)
%Overlays detected 2D segment obtained via subResSegment2DFit or any
%segment detection algorithm that complies with the subResSegment2D
%definition (see subResSegment2D.m for details)
%
%SYNPOSIS overlayFeaturesMovie(segmentParams,startend,saveMovie,movieName,...
%    filterSigma,showRaw,autoscaleImage)
%
%INPUT  segmentsInfo     : Cell array of size equal to the number of frames.
%                        Each cell contains the parameters of detection
%                        2D segments, i.e. an array of size Nx5, where N is
%                        the number of detected segments at a given frame
%                        and the number of columns corresponds to the xy
%                        position, amplitude, length and orientation of a
%                        segment.
%       startend       : Row vector indicating first and last frame to
%                        include in movie. Format: [startframe endframe].
%                        Optional. Default: [1 (maximum available frame)]
%       saveMovieInfo  : empty struct if no saving is requested. If saving is
%                        requested, structure with fields:
%           .dir          : Directory where the movie should be saved.
%           .filename     : Name of file where results should be saved.
%
%Sylvain Berlemont, 2010

%% input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--overSegment2DMovie: Incorrect number of input arguments!');
    return
end

%ask user for images
[fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image, even if not to be plotted');

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
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
    disp('--overSegment2DMovie: Bad file selection');
    return
    
end

%check startend and assign default if necessary
if nargin < 2 || isempty(startend)
    startend = [1 numFrames];
else
    startend(2) = min(startend(2),numFrames); %make sure that last frame does not exceed real last frame
end

%check whether to save movie
if nargin < 3 || isempty(saveMovieInfo)
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

%keep only the frames of interest
outFileList = outFileList(frame2fileMap(startend(1)):frame2fileMap(startend(2)));
frame2fileMap = frame2fileMap(startend(1):startend(2));
indxNotZero = find(frame2fileMap~=0);
frame2fileMap(indxNotZero) = frame2fileMap(indxNotZero) - frame2fileMap(indxNotZero(1)) + 1;

%initialize QT movie if it is to be saved
if saveMovie
    evalString = ['MakeQTMovie start ''' fullfile(saveMovieInfo.dir,saveMovieInfo.filename) ''''];
    eval(evalString);
end

%retain only the segmentParams of the frames of interest
if isempty(segmentParams)
    segmentParams = cell(startend(2)-startend(1)+1,1);
else
    segmentParams = segmentParams(startend(1):startend(2));
end

%get image size
imageRange = [1 isx; 1 isy];

%% make movie

%go over all specified frames
for iFrame = 1 : length(segmentParams)
    if frame2fileMap(iFrame) ~= 0 %if frame exists
        
        %read specified image
        imageStack = imread(outFileList{frame2fileMap(iFrame)});
        
    else %otherwise
        %make empty frame
        imageStack = zeros(isx,isy);
    end
    
    %plot image in current frame
    clf;
    
% axes('Position',[0 0 0.495 1]);
    imagesc(imageStack),colormap gray,axis image,axis off;
%     xlim(imageRange(2,:));
%     ylim(imageRange(1,:));
     hold on;
     textDeltaCoord = min(diff(imageRange,[],2))/20;
     text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
         textDeltaCoord,num2str(iFrame+startend(1)-1),'Color','white');
%     axes('Position',[0.505 0 0.495 1]);
%     imshow(imageStack,[]);
%     xlim(imageRange(2,:));
%     ylim(imageRange(1,:));
%     hold on;
    
    %plot segments
    if ~isempty(segmentParams{iFrame})
        xC = segmentParams{iFrame}(:,1);
        yC = segmentParams{iFrame}(:,2);
        l = segmentParams{iFrame}(:,4);
        t = segmentParams{iFrame}(:,5);

        line([xC - (l / 2) .* cos(t), xC + (l / 2) .* cos(t)]', ...
            [yC - (l / 2) .* sin(t), yC + (l / 2) .* sin(t)]', ...
            'Color', 'g');

        line(xC, yC, 'Color', 'g', 'Marker', '.', 'LineStyle', 'none');
    end
    
    %add frame to movie if movie is saved
    if saveMovie
        MakeQTMovie addfigure
    end
end

%finish movie
if saveMovie==1
    MakeQTMovie finish
end

%% ~~~ end ~~~