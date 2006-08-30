function [movieInfo,framesFailed,errFlag] = detectSubResFeatures2D_Movie(...
    movieParam,detectionParam)
%DETECTSUBRESFEATURES2D_MOVIE detects subresolution features in a series of images
%
%SYNOPSIS [detectedFeatures,errFlag] = detectSubResFeatures2D_Movie(movieParam,detectionParam)
%
%INPUT  movieParam    : Structure with fields
%           .imageDir     : Directory where images are stored
%           .candsDir     : Directory where cands (initial maxima) are stored.
%           .filenameBase : Filename base.
%           .firstImageNum: Numerical index of first image in movie.
%           .lastImageNum : Numerical index of last image in movie.
%           .digits4Enum  : Number of digits used to enumerate frames.
%       detectionParam: Structure with fields
%           .psfSigma     : Standard deviation of point spread function (in pixels).
%           .testAlpha    : Alpha-values for statistical tests. Optional.
%                           (See detectSubResFeatures2D for details). 
%           .visual       : 1 if user wants to view results; 0 otherwise. 
%                           Optional. Default: 0.
%           .doMMF        : 1 if user wants to do mixture-model fitting, 0
%                           otherwise. Optional. Default: 1.
%
%OUTPUT movieInfo     : Array of length "movie length" of structures 
%                       containing the fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       framesFailed  : Array indicating frames where detection has failed.
%                       If empty, then no frames failed.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, July 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movieInfo = [];
framesFailed = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--detectSubResFeatures2D_Movie: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get movie parameters
imageDir = movieParam.imageDir;
candsDir = movieParam.candsDir;
filenameBase = movieParam.filenameBase;
firstImageNum = movieParam.firstImageNum;
lastImageNum = movieParam.lastImageNum;

%get PSF sigma
psfSigma = detectionParam.psfSigma;

%get statistical test alpha values
if ~isfield(detectionParam,'testAlpha')
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05);
else
    testAlpha = detectionParam.testAlpha
end

%get visualization option
if ~isfield(detectionParam,'visual')
    visual = 0;
else
    visual = detectionParam.visual;
end

%check whether to do MMF
if ~isfield(detectionParam,'doMMF')
    doMMF = 1;
else
    doMMF = detectionParam.doMMF;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign leading zeros in numerical index of images
switch movieParam.digits4Enum
    case 4
        leadingZeros(1).value = '000';
        leadingZeros(2).value = '00';
        leadingZeros(3).value = '0';
        leadingZeros(4).value = '';
    case 3
        leadingZeros(1).value = '00';
        leadingZeros(2).value = '0';
        leadingZeros(3).value = '';
    case 2
        leadingZeros(1).value = '0';
        leadingZeros(2).value = '';
    case 1
        leadingZeros(1).value = '';
end

movieInfo(lastImageNum).xCoord = [];
movieInfo(lastImageNum).yCoord = [];
movieInfo(lastImageNum).amp = [];

%go through the images

for i=min(9999,lastImageNum):-1:max(1000,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(4).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(4).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,testAlpha,visual,doMMF);

        %save results
        movieInfo(i) = featuresInfo;
        
    catch %if detection fails
        
        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

end

for i=min(999,lastImageNum):-1:max(100,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(3).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(3).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,testAlpha,visual,doMMF);

        %save results
        movieInfo(i) = featuresInfo;

    catch %if detection fails

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

end

for i=min(99,lastImageNum):-1:max(10,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(2).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(2).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,testAlpha,visual,doMMF);

        %save results
        movieInfo(i) = featuresInfo;

    catch %if detection fails

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

end

for i=min(9,lastImageNum):-1:max(1,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(1).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(1).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,testAlpha,visual,doMMF);

        %save results
        movieInfo(i) = featuresInfo;

    catch %if detection fails

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

end


%%%%% ~~ the end ~~ %%%%%
