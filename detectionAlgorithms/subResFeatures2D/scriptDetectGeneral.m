
%% movie information
movieParam.imageDir = '/orchestra/groups/lccb-receptors/trialStuff/Goldstein/100406/'; %directory where images are
movieParam.filenameBase = 'APP_4_2_10_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 1; %number of last image in movie
movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.81; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 1; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 12; %Camera bit depth
detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

%% save results
saveResults.dir = '/orchestra/groups/lccb-receptors/trialStuff/Goldstein/100406/'; %directory where to save input and output
saveResults.filename = 'detectionTestAPP2.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
