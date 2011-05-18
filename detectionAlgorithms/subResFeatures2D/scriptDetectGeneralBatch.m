
%% define batch job locations

%image locations
imageDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/Martin/KJGermlings/G4/imagesCropped/',...
    };

%file name bases
filenameBase = {...
    'crop_G4_',...
    };

%directory for saving results
saveResDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/Martin/KJGermlings/G4/analysis/',...
    };

%background image locations
bgImageDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/Martin/KJGermlings/G4/bg/',...
    };

%background file name bases
bgFilenameBase = {...
    'crop_crop_G4_',...
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
        movieParam.firstImageNum = 1; %number of first image in movie
        movieParam.lastImageNum = 199; %number of last image in movie
        movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).
        
        %% detection parameters
        detectionParam.psfSigma = 1.6; %point spread function sigma (in pixels)
        detectionParam.testAlpha = struct('alphaR',1,'alphaA',1,'alphaD',1,'alphaF',0); %alpha-values for detection statistical tests
        detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
        detectionParam.doMMF = 0; %1 if mixture-model fitting, 0 otherwise
        detectionParam.bitDepth = 16; %Camera bit depth
        detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
        detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
        detectionParam.integWindow = 0; %number of frames before and after a frame for time integration
        
        %background info ...
        background.imageDir = bgImageDir{iMovie};
        background.filenameBase = bgFilenameBase{iMovie};
        background.alphaLocMaxAbs = 1e-15;
        detectionParam.background = background;
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionAll20.mat'; %name of file where input and output are saved
        % saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
