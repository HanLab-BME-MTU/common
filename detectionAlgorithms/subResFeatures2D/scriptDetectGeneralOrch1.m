
%% movie information
movieParam.imageDir = '/files/.retain-snapshots.d7d-w0d/LCCB/receptors/Galbraiths/data/lifeActAndCellEdge/120213/120213_Cs2C1_Lifeact/imagesLifeact/'; %directory where images are
movieParam.filenameBase = '120213_Cs2C1_CHO_mEos2LifeAct25ms_GFPfill5ms_'; %image file name base
movieParam.firstImageNum = 2; %number of first image in movie
movieParam.lastImageNum = 14400; %number of last image in movie
movieParam.digits4Enum = 5; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.15,'alphaD',0.15,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.15; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

detectionParam.calcMethod = 'g';

%absolute background info and parameters...
background.imageDir = '/files/.retain-snapshots.d7d-w0d/LCCB/receptors/Galbraiths/data/lifeActAndCellEdge/120213/120213_Cs2C1_Lifeact/bgLifeact/';
background.filenameBase = 'crop_120213_Cs2C1_CHO_mEos2LifeAct25ms_GFPfill5ms_';
background.alphaLocMaxAbs = 0.01;
detectionParam.background = background;

%% additional input

%save results
saveResults.dir = '/files/.retain-snapshots.d7d-w0d/LCCB/receptors/Galbraiths/data/lifeActAndCellEdge/120213/120213_Cs2C1_Lifeact/analysisLifeact/'; %directory where to save input and output
saveResults.filename = 'detectionAll1.mat'; %name of file where input and output are saved
% saveResults = 0;

%verbose state
verbose = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults,verbose);

%% exit

exit
