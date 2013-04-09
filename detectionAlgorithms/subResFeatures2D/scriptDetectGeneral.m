
%% movie information
movieParam.imageDir = 'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\121128_Cs1C3_Y773A\imagesAlphaVY773A\'; %directory where images are
movieParam.filenameBase = '121128_Cs1C3_CHO_mEos2AvBeta3Y773A_'; %image file name base
movieParam.firstImageNum = 2; %number of first image in movie
movieParam.lastImageNum = 9600; %number of last image in movie
movieParam.digits4Enum = 5; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

detectionParam.calcMethod = 'g';

%absolute background info and parameters...
background.imageDir = 'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\121128_Cs1C3_Y773A\bgAlphaVY773A\';
background.filenameBase = 'crop_121128_Cs1C3_CHO_mEos2AvBeta3Y773A_';
background.alphaLocMaxAbs = 0.01;
detectionParam.background = background;

%% additional input

%saveResults
saveResults.dir = 'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\121128_Cs1C3_Y773A\analysisAlphaVY773A\tracks\'; %directory where to save input and output
saveResults.filename = 'detectionAll1.mat'; %name of file where input and output are saved
% saveResults = 0;

%verbose state
verbose = 1;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults,verbose);
