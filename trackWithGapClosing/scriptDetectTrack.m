
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
movieParam.candsDir = '/mnt/sickkids/Fixed_Cy3movies_Nico/mov11Cy3_0030toodense&small/analysis_57to456/fsm/tack/cands/'; %directory where initial position estimates are
movieParam.imageDir = '/mnt/sickkids/Fixed_Cy3movies_Nico/mov11Cy3_0030toodense&small/cropped_57to456/'; %directory where images are
movieParam.filenameBase = 'crop_06-12-22_Calib_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 479; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%detectionParam
detectionParam.psfSigma = 0.511; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth

%saveResults
saveResults.dir = '/mnt/sickkids/Fixed_Cy3movies_Nico/mov11Cy3_0030toodense&small/analysis_57to456/mmf/'; %directory where to save input and output
saveResults.filename = 'detectionFixedCy3_mov11_003_57to456'; %name of file where input and output are saved

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
    'maxAmpRatio',10,...   %maximum ratio of amplitudes that allows the linking of 2 features
    'noLnkPrctl',-1);

%cost matrix for linking between frames again using statistical data on the tracks
costMatrices(2).costMatFun = 'costMatLogL_D2_V';
costMatrices(2).costMatParam = struct(...
    'cutoffProbD',0.9999,...    %cumulative probability of a square displacement beyond which linking is not allowed
    'cutoffProbA',1,...    %cumulative probability of an amplitude difference beyond which linking is not allowed
    'noLnkPrctl',-1,...
    'maxDist',3);               %search radius

%cost matrix for gap closing
costMatrices(3).costMatFun = 'costMatCloseGaps_D2_V';
costMatrices(3).costMatParam = struct(...
    'cutoffProbD1',0.9999,...   %cumulative probability of a square displacement beyond which gap closing is not allowed
    'cutoffProbA1',1,...   %cumulative probability of an amplitude difference beyond which gap closing is not allowed
    'cutoffProbD2',0.999,...    %cumulative probability of a square displacement beyond which merging/splitting are not allowed
    'cutoffProbA2',0.999,...    %cumulative probability of an amplitude difference beyond which merging/splitting are not allowed
    'noLnkPrctl',-1,...
    'maxDist',3*ones(20,1),...   %search radius for each time window
    'gapPenalty',zeros(20,1));   %penalty for closing gaps for each time window

%cost matrix for resolving merging and splitting conflicts
costMatrices(4).costMatFun = 'costVecLinkMS';
costMatrices(4).costMatParam = struct(...
    'cutoffProbD',0.9999,...    %cumulative probability of a square displacement beyond which merging/splitting are not allowed
    'cutoffProbA',0.9999);      %cumulative probability of an amplitude difference beyond which merging/splitting are not allowed

%gap closing parameters
gapCloseParam.timeWindow = 20;      %largest gap that can be closed
gapCloseParam.mergeSplit = 0;       %1 if merging/splitting are considered, 0 otherwise
gapCloseParam.segmentLength = 847;  %length of time segment for sequential gap closing

%iteration parameters
iterParam.tolerance = 0.05;   %maximum relative change of track statistical parameters to reach convergence
iterParam.lenFrac = 0.1;     %minimum length of tracks used for statistical analysis, as a function of the movie length

%saveResults
saveResults.dir = '/mnt/sickkids/Yoav_files/singleMoleculeYoav/calibration/analysis/mmf/'; %directory where to save input and output
saveResults.filename = 'trackedFeaturesCalibration'; %name of file where input and output are saved

%run the tracking function
[trackedFeatureNum,trackedFeatureInfo,errFlag] = trackWithGapClosing(...
    movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam,saveResults);

% %save the input and output
% save('trackedFeatures','costMatrices','gapCloseParam','iterParam',...
%     'trackedFeatureNum','trackedFeatureInfo');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

