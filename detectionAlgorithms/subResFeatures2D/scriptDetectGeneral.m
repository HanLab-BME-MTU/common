
%% movie information
movieParam.imageDir = '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20111018_LatranculinTreatment_TSP1_treated/Before_Lat/Images/'; %directory where images are
movieParam.filenameBase = 'Stream_unperturbed_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 100; %number of last image in movie
movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.001,'alphaA',0.4,'alphaD',1,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.25; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

%absolute background info and parameters...
background.imageDir = '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20111018_LatranculinTreatment_TSP1_treated/Before_Lat/bg/';
background.filenameBase = 'crop_Stream_unperturbed_';
background.alphaLocMaxAbs = 0.001;
detectionParam.background = background;

%% save results
saveResults.dir = '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20111018_LatranculinTreatment_TSP1_treated/Before_Lat/analysisKJ/'; %directory where to save input and output
saveResults.filename = 'detection11Frames001to100.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
