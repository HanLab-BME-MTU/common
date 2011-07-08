
%% movie information
movieParam.imageDir = '/home/kj35/orchestra/groups/lccb-receptors/codeTesting/'; %directory where images are
movieParam.filenameBase = 'test'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 1; %number of last image in movie
movieParam.digits4Enum = 1; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.5; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 1; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

%absolute background info and parameters...
% background.imageDir = '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/farnesylAndCellEdge/110519_CHO02/bgAlphaV01/';
% background.filenameBase = 'crop_110519_CHO02_mEos2Farn_';
% background.alphaLocMaxAbs = 0.01;
% detectionParam.background = background;

%% save results
% saveResults.dir = 'C:/Users/Anindya/Desktop/EAA1/Untreated/EAA/Analysis/'; %directory where to save input and output
% saveResults.filename = 'Untreated_EAA.mat'; %name of file where input and output are saved
saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
