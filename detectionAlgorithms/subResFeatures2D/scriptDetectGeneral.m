
%% movie information
movieParam.imageDir = '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_20/images/'; %directory where images are
movieParam.filenameBase = 'WT1_20_D3D_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 10; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 0.5; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.1,'alphaA',1,'alphaD',1,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 1; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

% %absolute background info and parameters...
background.imageDir = '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_20/bgImages/';
background.filenameBase = 'crop_WT1_20_D3D_';
background.alphaLocMaxAbs = 0.001;
detectionParam.background = background;

%% save results
saveResults.dir = '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_20/analysis/'; %directory where to save input and output
saveResults.filename = 'detectionTest1.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
