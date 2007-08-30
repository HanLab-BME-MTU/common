
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
%Get the sub-pixel positions of features - old
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Go to the directory ".../analysis/mmf"

%define the input variables

%movieParam
movieParam.candsDir = '/mnt/sickkids/Hiro/070817Cy3Control/C3con3/analysis/fsm/tack/cands/'; %directory where initial position estimates are
movieParam.imageDir = '/mnt/sickkids/Yoav/2007_08_20/monodisperse/100ms_cy3/cell_01/images/'; %directory where images are
movieParam.filenameBase = 'monodisperse_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 96; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%detectionParam
detectionParam.psfSigma = 1.7; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth

%saveResults
saveResults.dir = '/mnt/sickkids/Yoav/2007_08_20/monodisperse/100ms_cy3/cell_01/analysis/mmf/'; %directory where to save input and output
saveResults.filename = 'detectionNewFixedPSF2.mat'; %name of file where input and output are saved

%run the detection function
[movieInfo,emptyFrames,framesFailed] = detectSubResFeatures2D_Movie(...
    movieParam,detectionParam,saveResults);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the sub-pixel positions of features - new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Go to the directory ".../analysis/mmf"

%define the input variables

%movieParam
% movieParam.candsDir = '/mnt/sickkids/Hiro/070817Cy3Control/C3con3/analysis/fsm/tack/cands/'; %directory where initial position estimates are
movieParam.imageDir = '/mnt/sickkids/Yoav/2007_08_20/monodisperse/100ms_cy3/cell_01/images/'; %directory where images are
movieParam.filenameBase = 'monodisperse_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 96; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%detectionParam
detectionParam.psfSigma = 1.7; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth

%alphaLocMax
alphaLocMax = 0.005;

%saveResults
saveResults.dir = '/mnt/sickkids/Yoav/2007_08_20/monodisperse/100ms_cy3/cell_01/analysis/mmf/'; %directory where to save input and output
saveResults.filename = 'detectionNewFixedPSF2.mat'; %name of file where input and output are saved

%run the detection function
% [movieInfo,emptyFrames,framesFailed] = detectSubResFeatures2D_Movie(...
%     movieParam,detectionParam,saveResults);

[movieInfo,emptyFrames,framesFailedMMF,framesFailedLocMax] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,...
    alphaLocMax,saveResults);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track features over time - old
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stay in the directory ".../analysis/mmf"

%define the input variables

%cost matrix for the initial simple linking from frame to frame
costMatrices(1).costMatFun = 'costMatSimple';
costMatrices(1).costMatParam = struct('searchRadius',7,...   %maximum distance that allows the linking of 2 features
    'maxAmpRatio',10,...   %maximum ratio of amplitudes that allows the linking of 2 features
    'noLnkPrctl',-1);

%cost matrix for linking between frames again using statistical data on the tracks
costMatrices(2).costMatFun = 'costMatLogL_D2';
costMatrices(2).costMatParam = struct(...
    'cutoffProbD',0.9999,...    %cumulative probability of a square displacement beyond which linking is not allowed
    'cutoffProbA',1,...    %cumulative probability of an amplitude difference beyond which linking is not allowed
    'noLnkPrctl',-1,...
    'maxDist',7);               %search radius

%cost matrix for gap closing
costMatrices(3).costMatFun = 'costMatCloseGaps_D2';
costMatrices(3).costMatParam = struct(...
    'cutoffProbD1',0.9999,...   %cumulative probability of a square displacement beyond which gap closing is not allowed
    'cutoffProbA1',0.9999,...   %cumulative probability of an amplitude difference beyond which gap closing is not allowed
    'cutoffProbD2',0.999,...    %cumulative probability of a square displacement beyond which merging/splitting are not allowed
    'cutoffProbA2',0.999,...    %cumulative probability of an amplitude difference beyond which merging/splitting are not allowed
    'noLnkPrctl',-1,...
    'maxDist',7*ones(20,1),...   %search radius for each time window
    'gapPenalty',zeros(20,1));   %penalty for closing gaps for each time window

%cost matrix for resolving merging and splitting conflicts
costMatrices(4).costMatFun = 'costVecLinkMS';
costMatrices(4).costMatParam = struct(...
    'cutoffProbD',0.9999,...    %cumulative probability of a square displacement beyond which merging/splitting are not allowed
    'cutoffProbA',0.9999);      %cumulative probability of an amplitude difference beyond which merging/splitting are not allowed

%gap closing parameters
gapCloseParam.timeWindow = 20;      %largest gap that can be closed
gapCloseParam.mergeSplit = 0;       %1 if merging/splitting are considered, 0 otherwise
gapCloseParam.segmentLength = 687;  %length of time segment for sequential gap closing

%iteration parameters
iterParam.tolerance = 0.05;   %maximum relative change of track statistical parameters to reach convergence
iterParam.lenFrac = 0.1;     %minimum length of tracks used for statistical analysis, as a function of the movie length

%saveResults
saveResults.dir = '/mnt/sickkids/Yoav_files/2007_03_06_Calibration_slides/100ms_MaxExc_255Sens/analysis_57to456/mmf/'; %directory where to save input and output
saveResults.filename = 'tracksCalib_100ms_MaxExc_255Sens_57to456'; %name of file where input and output are saved

%run the tracking function
[trackedFeatureNum,trackedFeatureInfo,errFlag] = trackWithGapClosing(...
    movieInfo,costMatrices,'getTrackStats',gapCloseParam,iterParam,saveResults);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track features over time - new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stay in the directory ".../analysis/mmf"

%define the input variables

%some gap closing parameters
gapCloseParam.timeWindow = 10; %maximum allowed time gap (in frames) between a track end and a track start that allows linking them.
gapCloseParam.mergeSplit = 1; %1 if merging and splitting are considered, 0 if not.

%linking cost matrix parameters
%these are the parameters for linking detected features from one frame to
%the next in order to construct the initial tracks
costMatParam.minSearchRadiusL = 4; %minimum allowed search radius (in pixels). The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
costMatParam.maxSearchRadiusL = 4; %maximum allowed search radius (in pixels). Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
costMatParam.brownStdMultL = 3; %just keep this as 3. In In the final code I will probably hardwire this value in the code.
costMatParam.closestDistScaleL = 2; %same here. Keep as 2.
costMatParam.maxStdMultL = 20; %same here. Keep it as 20. I will explain these three parameters when I visit.

%gap closing cost matrix parameters
%these are the parameters for gap closing as well as merging and splitting
%these operations are performed on the initial tracks obtained in the
%previous step
costMatParam.minSearchRadiusCG = 4; %minimum allowed search radius (in pixels).
costMatParam.maxSearchRadiusCG = 4; %maximum allowed search radius (in pixels).
costMatParam.brownStdMultCG = 3*ones(gapCloseParam.timeWindow,1); %keep this as 3.
costMatParam.linStdMultCG = 3*ones(gapCloseParam.timeWindow,1); %keep this as 3.
costMatParam.timeReachConfB = 1; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed). I found a value of 1 works best, but you can play with this a little bit.
costMatParam.timeReachConfL = 10; %same as the previous parameter, but for the linear part of the motion. Again, I found that 5 works best, but you can play around with this parameter.
costMatParam.closestDistScaleCG = 2; %keep this as 2.
costMatParam.maxStdMultCG = 20; %and keep this as 20.
costMatParam.lenForClassify = 10; %keep this as 10.
costMatParam.maxAngle = 20; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
costMatParam.ampRatioLimitCG = [0.5000 2]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum the intensities of the 2 features that merge/split.

%parameters for using local density to expand search radius
%the search radius is generall calculated from the motion parameters of a
%feature or track. However, if a feature or track is really isolated from
%anything else, I expand its search radius.
useLocalDensity.link = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
useLocalDensity.cg = 1; %1 if you want to expand the search radius of isolated tracks in the gap closing step.
useLocalDensity.nnWindowL = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
useLocalDensity.nnWindowCG = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look to see a track's nearest neighbor at its end/start (in the gap closing step).

%saveResults
saveResults.dir = '/mnt/sickkids/Hiro/070820Cy3Titration/80000_30_6/analysis/mmf'; %directory where to save input and output
saveResults.filename = 'tracks_070820Cy3Titration_80000_30_6.mat'; %name of file where input and output are saved

%run the tracking function
[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalman(movieInfo,...
    costMatParam,gapCloseParam,[],useLocalDensity,saveResults);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

