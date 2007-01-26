
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the necessary directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%In the movie directory of interest, make 3 directories: "images",
%"cropped" and "analysis". 

%In the directory "analysis", make 2 directories: "fsm" and "mmf".


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the initial estimates of feature positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Go to the directory ".../analysis/fsm" and type fsmCenter in command line

%Put the proper movie parameters

%Set up project

%crop a background part of the image and estimate noise parameters

%run speckTackle to get initial feature positions as "cands" files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the sub-pixel positions of features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Go to the directory ".../analysis/mmf"

%define the input variables

%movieParam
movieParam.candsDir = '/mnt/mit/khuloudDir/synapse/July2006/0707_ML7/dish3/ML7_1-good/analysisPart/fsm/tack/cands/'; %directory where initial position estimates are
movieParam.imageDir = '/mnt/mit/khuloudDir/synapse/July2006/0707_ML7/dish3/ML7_1-good/imagesPart/'; %directory where images are
movieParam.filenameBase = 'crop_070706_dish3_ML7_105frps_gain205_1_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 1; %number of last image in movie
movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

%detectionParam
detectionParam.psfSigma = 1.078; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 14; %Camera bit depth

%saveResults
saveResults.dir = '/mnt/mit/khuloudDir/synapse/July2006/0707_ML7/dish3/ML7_1-good/analysisPart/mmf/'; %directory where to save input and output
saveResults.filename = 'detectedFeatures'; %name of file where input and output are saved

%run the detection function
[movieInfo,emptyFrames,framesFailed] = detectSubResFeatures2D_Movie(...
    movieParam,detectionParam,saveResults);

% %save the input and output
% save('detectedFeatures','movieParam','detectionParam','movieInfo','emptyFrames','framesFailed');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track features over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stay in the directory ".../analysis/mmf"

%define the input variables

%cost matrix for the initial simple linking from frame to frame
costMatrices(1).costMatFun = 'costMatSimple';
costMatrices(1).costMatParam = struct('searchRadius',3,...   %maximum distance that allows the linking of 2 features
    'maxAmpRatio',2,...   %maximum ratio of amplitudes that allows the linking of 2 features
    'noLnkPrctl',-1);

%cost matrix for linking between frames again using statistical data on the tracks
costMatrices(2).costMatFun = 'costMatLogL';
costMatrices(2).costMatParam = struct(...
    'cutoffProbD',0.9999,...    %cumulative probability of a square displacement beyond which linking is not allowed
    'cutoffProbA',0.9999,...    %cumulative probability of an amplitude difference beyond which linking is not allowed
    'noLnkPrctl',-1);

%cost matrix for gap closing
costMatrices(3).costMatFun = 'costMatCloseGaps';
costMatrices(3).costMatParam = struct(...
    'cutoffProbD1',0.9999,...   %cumulative probability of a square displacement beyond which gap closing is not allowed
    'cutoffProbA1',0.9999,...   %cumulative probability of an amplitude difference beyond which gap closing is not allowed
    'cutoffProbD2',0.999,...   %cumulative probability of a square displacement beyond which merging/splitting are not allowed
    'cutoffProbA2',0.999,...   %cumulative probability of an amplitude difference beyond which merging/splitting are not allowed
    'noLnkPrctl',-1);

%cost matrix for resolving merging and splitting conflicts
costMatrices(4).costMatFun = 'costVecLinkMS';
costMatrices(4).costMatParam = struct(...
    'cutoffProbD',0.9999,...    %cumulative probability of a square displacement beyond which merging/splitting are not allowed
    'cutoffProbA',0.9999);      %cumulative probability of an amplitude difference beyond which merging/splitting are not allowed

%gap closing parameters
gapCloseParam.timeWindow = 5;       %largest gap that can be closed
gapCloseParam.mergeSplit = 0;       %1 if merging/splitting are considered, 0 otherwise
gapCloseParam.segmentLength = 50;  %length of time segment for sequential gap closing

%iteration parameters
iterParam.tolerance = 0.05;   %maximum relative change of track statistical parameters to reach convergence
iterParam.lenFrac = 0.5;      %minimum length of tracks used for statistical analysis, as a function of the movie length

%saveResults
saveResults.dir = '/mnt/mit/khuloudDir/synapse/July2006/0707_ML7/dish3/ML7_1-good/analysisPart/mmf/'; %directory where to save input and output
saveResults.filename = 'trackedFeatures'; %name of file where input and output are saved

%run the tracking function
[trackedFeatureNum,trackedFeatureInfo,errFlag] = trackWithGapClosing(...
    movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam,saveResults);

% %save the input and output
% save('trackedFeatures','costMatrices','gapCloseParam','iterParam',...
%     'trackedFeatureNum','trackedFeatureInfo');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

