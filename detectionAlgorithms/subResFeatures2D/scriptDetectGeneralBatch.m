
%% define batch job locations

%image locations
imageDir = {...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131127\imagesAlphaVY773A\',...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131202\imagesAlphaVY773A\',...
    };

%file name bases
filenameBase = {...
    '131127_Cs3C2_CHO_Y773A_',...
    '131202_Cs3C1_CHO_Y773A_',...
    };

%directory for saving results
saveResDir = {...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131127\analysisAlphaVY773A\',...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131202\analysisAlphaVY773A\',...
    };

%background image locations
bgImageDir = {...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131127\bgAlphaVY773A\',...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131202\bgAlphaVY773A\',...
    };

%background file name bases
bgFilenameBase = {...
    'crop_131127_Cs3C2_CHO_Y773A_',...
    'crop_131202_Cs3C1_CHO_Y773A_',...
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
        movieParam.lastImageNum = 14400; %number of last image in movie
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
        
        %background info ...
        background.imageDir = bgImageDir{iMovie};
        background.filenameBase = bgFilenameBase{iMovie};
        background.alphaLocMaxAbs = 0.01;
        detectionParam.background = background;
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionAll1.mat'; %name of file where input and output are saved
        %         saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
