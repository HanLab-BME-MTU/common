function [movieInfo,emptyFrames,framesFailedMMF,framesFailedLocMax,...
    psfSigma] = detectSubResFeatures2D_StandAlone(movieParam,...
    detectionParam,alphaLocMax,estimateSigma,saveResults)
%DETECTSUBRESFEATURES2D_STANDALONE detects subresolution features in a series of images
%
%SYNOPSIS [movieInfo,emptyFrames,framesFailedMMF,framesFailedLocMax,...
%    psfSigma] = detectSubResFeatures2D_StandAlone(movieParam,...
%    detectionParam,alphaLocMax,saveResults)
%
%INPUT  movieParam    : Structure with fields
%           .imageDir     : Directory where images are stored
%           .filenameBase : Filename base.
%           .firstImageNum: Numerical index of first image in movie.
%           .lastImageNum : Numerical index of last image in movie.
%           .digits4Enum  : Number of digits used to enumerate frames.
%       detectionParam: Structure with fields
%           .psfSigma     : Initial guess for standard deviation of point
%                           spread function (in pixels).
%           .testAlpha    : Alpha-values for statistical tests in 
%                           detectSubResFeatures2D. Optional.
%                           (See detectSubResFeatures2D for details).
%           .visual       : 1 if user wants to view results; 0 otherwise.
%                           Optional. Default: 0.
%           .doMMF        : 1 if user wants to do mixture-model fitting, 0
%                           otherwise. Optional. Default: 1.
%           .bitDepth     : Camera bit depth. Optional. Default: 14.
%       alphaLocMax   : Alpha value for statistical test in local maxima
%                       detection. Optional. default: 0.005.
%       estimateSigma : 1 to refine PSF sigma value using the data, 0
%                       otherwise. Optional. Default: 1.
%       saveResults   : 0 if no saving is requested. 
%                       If saving is requested, structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                           Optional. Default: detectedFeatures.
%                       Whole structure optional.
%
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT movieInfo     : Array of length "movie length" of structures
%                       containing the fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       emptyFrames   : Array indicating frames where no features were
%                       detected.
%       framesFailedMMF: Array indicating frames where mixture-model fitting
%                        has failed.
%       framesFailedLocMax: Array indicating frames where initial detection
%                           of local maxima has failed.
%       psfSigma      : Standard deviation of point spread function as
%                       estimated from the data.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, July 2006

%% Output

movieInfo = [];
emptyFrames = [];
framesFailedMMF = [];
framesFailedLocMax = [];
psfSigma = [];

%% Input + Pre-processing

%check whether correct number of input arguments was used
if nargin < 2
    disp('--detectSubResFeatures2D_StandAlone: Incorrect number of input arguments!');
    return
end

%get movie parameters
imageDir = movieParam.imageDir;
filenameBase = movieParam.filenameBase;
firstImageNum = movieParam.firstImageNum;
lastImageNum = movieParam.lastImageNum;
digits4Enum = movieParam.digits4Enum;

%get initial guess of PSF sigma
psfSigma = detectionParam.psfSigma;

%get statistical test alpha values
if ~isfield(detectionParam,'testAlpha') || isempty(detectionParam.testAlpha)
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0.05);
else
    testAlpha = detectionParam.testAlpha;
end

%get visualization option
if ~isfield(detectionParam,'visual') || isempty(detectionParam.visual)
    visual = 0;
else
    visual = detectionParam.visual;
end

%check whether to do MMF
if ~isfield(detectionParam,'doMMF') || isempty(detectionParam.doMMF)
    doMMF = 1;
else
    doMMF = detectionParam.doMMF;
end

%get camera bit depth
if ~isfield(detectionParam,'bitDepth') || isempty(detectionParam.bitDepth)
    bitDepth = 14;
else
    bitDepth = detectionParam.bitDepth;
end

%get alpha-value for local maxima detection
if nargin < 3 || isempty(alphaLocMax)
    alphaLocMax = 0.005;
end

%check whether to estimate PSF sigma from the data
if nargin < 4 || isempty(estimateSigma)
    estimateSigma = 1;
end
    
%determine where to save results
if nargin < 5 || isempty(saveResults) %if nothing was input
    saveResDir = pwd;
    saveResFile = 'detectedFeatures';
    saveResults.dir = pwd;
else
    if isstruct(saveResults)
        if ~isfield(saveResults,'dir') || isempty(saveResults.dir)
            saveResDir = pwd;
        else
            saveResDir = saveResults.dir;
        end
        if ~isfield(saveResults,'filename') || isempty(saveResults.filename)
            saveResFile = 'detectedFeatures';
        else
            saveResFile = saveResults.filename;
        end
    else
        saveResults = 0;
    end
end

%store the string version of the numerical index of each image
enumString = getStringIndx(digits4Enum);

%% Image reading + filtering

%turn warnings off
warningState = warning('off','all');

%get image related parameters
imageTmp = imread([imageDir filenameBase enumString(1,:) '.tif']); %first image
[imageSizeX,imageSizeY] = size(imageTmp); %image size
imageIndx = firstImageNum : lastImageNum; %image indices
numImages = lastImageNum - firstImageNum + 1; %number of images
clear imageTmp

%initialize progress display
progressText(0,'Reading images');

%read images 
image = zeros(imageSizeX,imageSizeY,numImages);
for iImage = 1 : numImages
    
    %store images in array
    image(:,:,iImage) = imread([imageDir filenameBase enumString(imageIndx(iImage),:) '.tif']);    

    %display progress
    progressText(iImage/numImages,'Reading images');

end

%normalize images
image = image / (2^bitDepth-1);

%initialize progress display
progressText(0,'Filtering images');

%filter images
imageF = zeros(imageSizeX,imageSizeY,numImages);
for iImage = 1 : numImages

    imageF(:,:,iImage) = Gauss2D(image(:,:,iImage),psfSigma);

    %display progress
    progressText(iImage/numImages,'Filtering images');
end

%% Background noise estimation

%get the intensities in the last 5 images
last5start = max(numImages-5,1);
imageLast5 = double(image(:,:,last5start:numImages));

%estimate the background noise mean and standard deviation
%use robustMean to get mean and std of intensities
%in this method, the intensities of actual features will look like
%outliers, so we are effectively getting the mean and std of the background
[bgMean,bgStd] = robustMean(imageLast5(:));

%do the same for the filtered image
imageLast5 = double(imageF(:,:,last5start:numImages));
[bgMeanF,bgStdF] = robustMean(imageLast5(:));


%% Local maxima detection

%initialize structure saving local maxima information
localMaxima = repmat(struct('cands',[]),numImages,1);

%initialize progress display
progressText(0,'Detecting local maxima');

%go over all images ...
for iImage = 1 : numImages

    try

        %call locmax2d to get local maxima in filtered image
        fImg = locmax2d(imageF(:,:,iImage),[1 1]*ceil(3*psfSigma));
        
        %get positions and amplitudes of local maxima
        [localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);

        %calculate the p-value corresponding to the local maxima's amplitudes
        %assume that background intensity in filtered image is normally
        %distributed with mean bgMeanF and standard deviation bgStdF
        pValue = 1 - normcdf(localMaxAmp,bgMeanF,bgStdF);

        %retain only those maxima with significant amplitude
        keepMax = find(pValue < alphaLocMax);
        localMaxPosX = localMaxPosX(keepMax);
        localMaxPosY = localMaxPosY(keepMax);
        localMaxAmp = localMaxAmp(keepMax);
        pValue = pValue(keepMax);
        numLocalMax = length(keepMax);

        %construct cands structure
        if numLocalMax == 0 %if there are no local maxima

            cands = [];
            emptyFrames = [emptyFrames; imageIndx(iImage)];

        else %if there are local maxima

            %define background mean and status
            cands = repmat(struct('IBkg',bgMean,'status',1,...
                'Lmax',[],'amp',[],'pValue',[]),numLocalMax,1);
            
            %store maxima positions, amplitudes and p-values
            for iMax = 1 : numLocalMax
                cands(iMax).Lmax = [localMaxPosX(iMax) localMaxPosY(iMax)];
                cands(iMax).amp = localMaxAmp(iMax);
                cands(iMax).pValue = pValue(iMax);
            end

        end

        %add the cands of the current image to the rest
        localMaxima(iImage).cands = cands;
         
    catch

        %if local maxima detection fails, make cands empty
        localMaxima(iImage).cands = [];
        
        %add this frame to the array of frames with failed local maxima
        %detection and to the array of empty frames
        framesFailedLocMax = [framesFailedLocMax; imageIndx(iImage)];
        emptyFrames = [emptyFrames; imageIndx(iImage)];
        
    end

    %display progress
    progressText(iImage/numImages,'Detecting local maxima');

end

%make a list of images that have local maxima
goodImages = setxor(1:numImages,emptyFrames);

%% PSF sigma estimation

if estimateSigma

    %specify which parameters to fit for sigma estimation
    fitParameters = [{'X1'} {'X2'} {'A'} {'Sxy'} {'B'}];

    %save input PSF sigma in new variable and empty psfSigma for estimation
    psfSigma0 = psfSigma;
    psfSigma = [];
    
    %calculate some numbers that get repeated many times
    psfSigma5 = ceil(5*psfSigma0);

    %initialize progress display
    progressText(0,'Estimating PSF sigma');

    %go over the first 10 images and find isolated features
%     for iImage =  1 : min(numImages,10)
    for iImage =  1 : numImages

        %get feature positions and amplitudes
        featPos = vertcat(localMaxima(iImage).cands.Lmax);
        featAmp = vertcat(localMaxima(iImage).cands.amp);

        %retain only features that are more than 5*psfSigma0 away from boundaries
        feat2use = find(featPos(:,1) > psfSigma5 & ...
            featPos(:,1) < imageSizeX - psfSigma5 & ...
            featPos(:,2) > psfSigma5 & featPos(:,2) < imageSizeY - psfSigma5);
        featPos = featPos(feat2use,:);
        featAmp = featAmp(feat2use);

        %find nearest neighbor distances
        nnDist = createDistanceMatrix(featPos,featPos);
        nnDist = sort(nnDist,2);
        nnDist = nnDist(:,2);

        %retain only features whose nearest neighbor is more than 20*psfSigma0
        %away
        feat2use = find(nnDist > ceil(10*psfSigma0));
        featPos = featPos(feat2use,:);
        featAmp = featAmp(feat2use);

        %retain only features with amplitudes between the 25th and 75th
        %percentiles
        percentile25 = prctile(featAmp,25);
        percentile75 = prctile(featAmp,75);
        feat2use = find(featAmp > percentile25 & featAmp < percentile75);
        featPos = featPos(feat2use,:);
        featAmp = featAmp(feat2use);

        %go over the selected features and estimate psfSigma
        numFeats = length(featAmp);
        parameters = zeros(numFeats,5);
        for iFeat = 1 : numFeats

            %crop image around selected feature
            lowerBound = featPos(iFeat,:) - psfSigma5;
            upperBound = featPos(iFeat,:) + psfSigma5;
            imageCropped = image(lowerBound(1):upperBound(1),...
                lowerBound(2):upperBound(2),iImage);

            %make initial guess for fit (in the order given in fitParameters)
            initGuess = [psfSigma5+1 psfSigma5+1 featAmp(iFeat) psfSigma0 bgMean];

            %fit image and estimate sigma of Gaussian
            parameters(iFeat,:) = GaussFitND(imageCropped,[],fitParameters,initGuess);

        end

        %add to array of sigmas
        psfSigma = [psfSigma; parameters(:,4)];

    %display progress
%     progressText(iImage/min(numImages,10),'Estimating PSF sigma');
    progressText(iImage/numImages,'Estimating PSF sigma');

    end

    %esimate psfSigma as the robust mean of all the sigmas from the fits
    numCalcs = length(psfSigma);
    [psfSigma,sigmaStd,inlierIndx] = robustMean(psfSigma);
    
    disp(sprintf('PSF sigma = %1.3f (%d inliers out of %d observations)',...
        psfSigma,length(inlierIndx),numCalcs));
    
end

%% Mixture-model fitting

%initialize movieInfo
clear movieInfo
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numImages,1);

%initialize progress display
progressText(0,'Mixture-model fitting');

%go over all non-empty images ...
for iImage = goodImages

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image(:,:,iImage),...
            localMaxima(iImage).cands,psfSigma,testAlpha,visual,doMMF,1,0,bgStd);

        %save results
        movieInfo(iImage) = featuresInfo;

        %check whether frame is empty
        if isempty(featuresInfo.xCoord)
            emptyFrames = [emptyFrames; imageIndx(iImage)];
        end

    catch %if detection fails

        %label frame as empty
        emptyFrames = [emptyFrames; imageIndx(iImage)];

        %add this frame to the array of frames with failed mixture-model
        %fitting
        framesFailedMMF = [framesFailedMMF; imageIndx(iImage)];
        
    end

    %display progress
    progressText(iImage/numImages,'Mixture-model fitting');

end

%% Post-processing

%sort list of empty frames
emptyFrames = sort(emptyFrames);

%indicate correct frames in movieInfo
tmptmp = movieInfo;
clear movieInfo
movieInfo(firstImageNum:lastImageNum) = tmptmp;

%save results
if isstruct(saveResults)
    save([saveResDir filesep saveResFile],'movieParam','detectionParam',...
        'alphaLocMax','estimateSigma','movieInfo','emptyFrames',...
        'framesFailedMMF','framesFailedLocMax','psfSigma');
end

%go back to original warnings state
warning(warningState);

%%


%% Subfunctions

function enumString = getStringIndx(digits4Enum)

switch digits4Enum
    case 4
        enumString = repmat('0',9999,4);
        for i = 1 : 9
            enumString(i,:) = ['000' num2str(i)];
        end
        for i = 10 : 99
            enumString(i,:) = ['00' num2str(i)];
        end
        for i = 100 : 999
            enumString(i,:) = ['0' num2str(i)];
        end
        for i = 1000 : 9999
            enumString(i,:) = num2str(i);
        end
    case 3
        enumString = repmat('0',999,3);
        for i = 1 : 9
            enumString(i,:) = ['00' num2str(i)];
        end
        for i = 10 : 99
            enumString(i,:) = ['0' num2str(i)];
        end
        for i = 100 : 999
            enumString(i,:) = num2str(i);
        end
    case 2
        enumString = repmat('0',99,2);
        for i = 1 : 9
            enumString(i,:) = ['0' num2str(i)];
        end
        for i = 10 : 99
            enumString(i,:) = num2str(i);
        end
    case 1
        enumString = repmat('0',9,1);
        for i = 1 : 9
            enumString(i,:) = num2str(i);
        end
end
