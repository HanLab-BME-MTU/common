
%% movie information
movieParam.imageDir = '/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/100323_Cs1_CHO05/imagesAlphaV/'; %directory where images are
movieParam.filenameBase = '100323_Cs1_CHO05m_mEosAv_10msEGFPfill_7p2K_25ms_'; %image file name base
movieParam.firstImageNum = 100; %number of first image in movie
movieParam.lastImageNum = 7200; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.01,'alphaA',0.01,'alphaD',0.01,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

%% save results
saveResults.dir = '/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/100323_Cs1_CHO05/analysisAlphaV/'; %directory where to save input and output
saveResults.filename = 'detection2Frames0100to7200.mat'; %name of file where input and output are saved
% saveResults = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
