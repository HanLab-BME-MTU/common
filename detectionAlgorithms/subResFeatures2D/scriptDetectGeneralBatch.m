
%% define batch job locations

%image locations
imageDir = {...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C1_CHO_Farn/imagesFarn01/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C2_CHO_Farn/imagesFarn01/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C4_CHO_Farn/imagesFarn01/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs2C1_CHO_Farn/imagesFarn01/',...
    };

%file name bases
filenameBase = {...
    '110829_Cs1C1_CHO_mEos2Farn_',...
    '110829_Cs1C2_CHO_mEos2Farn_',...
    '110829_Cs1C4_CHO_mEos2Farn_',...
    '110829_Cs2C1_CHO_mEos2Farn_',...
    };

%directory for saving results
saveResDir = {...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C1_CHO_Farn/analysisFarn/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C2_CHO_Farn/analysisFarn/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C4_CHO_Farn/analysisFarn/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs2C1_CHO_Farn/analysisFarn/',...
    };

%background image locations
bgImageDir = {...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C1_CHO_Farn/bg01/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C2_CHO_Farn/bg01/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C4_CHO_Farn/bg01/',...
    '/home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs2C1_CHO_Farn/bg01/',...
    };

%background file name bases
bgFilenameBase = {...
    'crop_110829_Cs1C1_CHO_mEos2Farn_',...
    'crop_110829_Cs1C2_CHO_mEos2Farn_',...
    'crop_110829_Cs1C4_CHO_mEos2Farn_',...
    'crop_110829_Cs2C1_CHO_mEos2Farn_',...
    };

%% calculate number of movies
numMovies = length(filenameBase);

for iMovie = 1 : numMovies
    
    try
        
        %display movie number
        disp(['Movie ' num2str(iMovie) ' / ' num2str(numMovies) ' ...'])
        
        %% movie information
        movieParam.imageDir = imageDir{iMovie}; %directory where images are
        movieParam.filenameBase = filenameBase{iMovie}; %image file name base
        movieParam.firstImageNum = 2; %number of first image in movie
        movieParam.lastImageNum = 2400; %number of last image in movie
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
        
        %background info ...
        background.imageDir = bgImageDir{iMovie};
        background.filenameBase = bgFilenameBase{iMovie};
        background.alphaLocMaxAbs = 0.01;
        detectionParam.background = background;
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionAll01.mat'; %name of file where input and output are saved
        %         saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
