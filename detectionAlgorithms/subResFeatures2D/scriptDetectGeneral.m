
%% movie information
movieParam.imageDir = '/home/kj35/orchestra/groups/lccb-receptors/Grinstein/Nico/1011_MTsAndCD36/images/'; %directory where images are
movieParam.filenameBase = 'movie_13_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 500; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 2.1; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.01; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 1; %number of frames before and after a frame for time integration

% %absolute background info and parameters...
% background.imageDir = '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_20/bgImages/';
% background.filenameBase = 'crop_WT1_20_D3D_';
% background.alphaLocMaxAbs = 0.001;
% detectionParam.background = background;

%% save results
saveResults.dir = '/home/kj35/orchestra/groups/lccb-receptors/Grinstein/Nico/1011_MTsAndCD36/analysis/'; %directory where to save input and output
saveResults.filename = 'detectionAll8.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
