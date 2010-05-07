
%% movie information
movieParam.imageDir = '???\CD36\con_cy3_02\images4detection\'; %directory where images are
movieParam.filenameBase = 'crop_071017_37CLNB_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 20; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% save results
saveResults.dir = '???\CD36\con_cy3_02\analysis\'; %directory where to save input and output
saveResults.filename = 'detection20Frames_1.mat'; %name of file where input and output are saved

%% detection parameters
detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
detectionParam.testAlpha = struct('alphaR',0.01,'alphaA',0.01,'alphaD',0.01,'alphaF',0); %alpha-values for detection statistical tests

detectionParam.psfSigma = 2.1; %point spread function sigma (in pixels)
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 1; %number of frames before and after a frame for time integration

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
