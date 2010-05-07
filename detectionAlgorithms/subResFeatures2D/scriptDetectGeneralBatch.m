
%% define batch job locations

%image locations
imageDir = {...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_13/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_26/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_34/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_41/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_42/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_8/images/',...
    };

%file name bases
filenameBase = {...
    'img_',...
    'img_',...
    'img_',...
    'img_',...
    'img_',...
    'img_',...
    };

%directory for saving results
saveResDir = {...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_13/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_26/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_34/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_41/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_42/images/',...
    '/orchestra/groups/lccb-receptors/Grinstein/Hiro/20100429_Qdot655BrownianLinear/export/1_200_8/images/',...
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
        movieParam.lastImageNum = 315; %number of last image in movie
        movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).
        
        %% detection parameters
        detectionParam.psfSigma = 2.1; %point spread function sigma (in pixels)
        detectionParam.testAlpha = struct('alphaR',0.00001,'alphaA',0.00001,'alphaD',0.00001,'alphaF',0); %alpha-values for detection statistical tests
        detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
        detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
        detectionParam.bitDepth = 16; %Camera bit depth
        detectionParam.alphaLocMax = 0.001; %alpha-value for initial detection of local maxima
        detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
        detectionParam.integWindow = 0; %number of frames before and after a frame for time integration
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionAll4.mat'; %name of file where input and output are saved
        % saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
