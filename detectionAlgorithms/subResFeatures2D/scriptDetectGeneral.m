
%% movie information
movieParam.imageDir = '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_06/imagesCropped/'; %directory where images are
movieParam.filenameBase = 'crop_071017_37CLNB_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 99; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 2.0; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 1; %number of frames before and after a frame for time integration

%% save results
saveResults.dir = '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_06/analysisCropped/'; %directory where to save input and output
saveResults.filename = 'detection_100216_afterBugFix_2.mat'; %name of file where input and output are saved
% saveResults = 1;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
