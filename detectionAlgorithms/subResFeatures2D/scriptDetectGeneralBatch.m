
%% define batch job locations

%image locations
imageDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03A/imagesAlphaV/',...
    '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03B/imagesAlphaV/',...
    };

%file name bases
filenameBase = {...
    '110114_Cs1_CHO03A_',...
    '110114_Cs1_CHO03B_',...
    };

%directory for saving results
saveResDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03A/analysisAlpha/',...
    '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03B/analysisAlpha/',...
    };

%background image locations
bgImageDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03A/bgAlphaV/',...
    '/home/kj35/orchestra/groups/lccb-receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03B/bgAlphaV/',...
    };

%background file name bases
bgFilenameBase = {...
    'crop_110114_Cs1_CHO03A_',...
    'crop_110114_Cs1_CHO03B_',...
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
        movieParam.lastImageNum = 6800; %number of last image in movie
        movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).
        
        %% detection parameters
        detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
        detectionParam.testAlpha = struct('alphaR',0.001,'alphaA',0.01,'alphaD',0.001,'alphaF',0); %alpha-values for detection statistical tests
        detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
        detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
        detectionParam.bitDepth = 16; %Camera bit depth
        detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
        detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
        detectionParam.integWindow = 0; %number of frames before and after a frame for time integration
        
        %background info ...
        background.imageDir = bgImageDir{iMovie};
        background.filenameBase = bgFilenameBase{iMovie};
        background.alphaLocMaxAbs = 0.01;
        detectionParam.background = background;
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionAll2.mat'; %name of file where input and output are saved
        % saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
