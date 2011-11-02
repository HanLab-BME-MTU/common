
%% movie information
movieParam.imageDir = '/home/kj35/files/LCCB/receptors/codeTesting/Olivo-Marin/trackingPerformanceEvaluation/real/tiffFiles/'; %directory where images are
movieParam.filenameBase = 'im_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 128; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.85; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.01,'alphaA',0.25,'alphaD',0.01,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = [0.25 0.25]; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = [0 1]; %number of frames before and after a frame for time integration

% %absolute background info and parameters...
% background.imageDir = '/home/kj35/files/LCCB/maki/newUnarchived/111101_EdwardHarry_1sKinetochoreMovies/bg3/';
% background.filenameBase = 'projection_3_';
% background.alphaLocMaxAbs = 0.01;
% detectionParam.background = background;

%% save results
saveResults.dir = '/home/kj35/files/LCCB/receptors/codeTesting/Olivo-Marin/trackingPerformanceEvaluation/real/analysis/'; %directory where to save input and output
saveResults.filename = 'detectionAll9.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
