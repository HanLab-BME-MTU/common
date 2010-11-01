
%% movie information
movieParam.imageDir = '/home/kj35/orchestra/groups/lccb-receptors/Grinstein/Hiro/080421Qdot_8ms/Q2/images/'; %directory where images are
movieParam.filenameBase = '0407_C3Qd_Jasp_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 200; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 2.1; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.001; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

% %absolute background info and parameters...
% background.imageDir = '/home/kj35/orchestra/groups/lccb-receptors/Grinstein/Hiro/101021_OldTitrationData/0p001/imagesRenamed/';
% background.filenameBase = 'crop_2010_07_28_50_';
% background.alphaLocMaxAbs = 0.001;
% detectionParam.background = background;

%% save results
saveResults.dir = '/home/kj35/orchestra/groups/lccb-receptors/Grinstein/Hiro/080421Qdot_8ms/Q2/analysis/'; %directory where to save input and output
saveResults.filename = 'detectionFrame1to100_101021.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
