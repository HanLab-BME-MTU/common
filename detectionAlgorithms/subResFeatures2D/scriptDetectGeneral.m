
%% movie information
movieParam.imageDir = 'U:\Hiro\090310_FixedFabDilution\0p3\'; %directory where images are
movieParam.filenameBase = '0p3_100_'; %image file name base
movieParam.firstImageNum = 2; %number of first image in movie
movieParam.lastImageNum = 2; %number of last image in movie
movieParam.digits4Enum = 1; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.7; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 1; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

%% save results
% saveResults.dir = '/mnt/sickkids/codeTesting/testDetection/GapLengthDistributionCSHL/'; %directory where to save input and output
% saveResults.filename = 'detectionAlpha0p075.mat'; %name of file where input and output are saved
saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
