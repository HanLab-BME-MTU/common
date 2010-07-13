
%% movie information
movieParam.imageDir = '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/100618/CHO09/imagesAlphaV/'; %directory where images are
movieParam.filenameBase = '100618_Cs3_CHO09_mEosPax_10GFPfill_8202_25_MOVED_'; %image file name base
movieParam.firstImageNum = 2; %number of first image in movie
movieParam.lastImageNum = 8200; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.001,'alphaA',0.01,'alphaD',0.001,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = [0.05 0.1]; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = [0 1]; %number of frames before and after a frame for time integration

%% save results
saveResults.dir = '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/100618/CHO09/analysisAlphaV/'; %directory where to save input and output
saveResults.filename = 'detectionAll2.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
