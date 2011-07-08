
%% define batch job locations

%image locations
imageDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/endoCells/AnindyaChanda/20110702_Colocalization/Unstimulated/Field4/',...
    };

%file name bases
filenameBase = {...
    'experiment7_w3642_',...
    };

%directory for saving results
saveResDir = {...
    '/home/kj35/orchestra/groups/lccb-receptors/endoCells/AnindyaChanda/20110702_Colocalization/Unstimulated/Field4/',...
    };

% %background image locations
% bgImageDir = {...
%     '/home/kj35/orchestra/groups/lccb-receptors/Martin/KJWildTypeHyphae/WT1/bg/',...
%     };
% 
% %background file name bases
% bgFilenameBase = {...
%     'crop_WT1_',...
%     };

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
        movieParam.lastImageNum = 1; %number of last image in movie
        movieParam.digits4Enum = 2; %number of digits used for frame enumeration (1-4).
        
        %% detection parameters
        detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
        detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.1,'alphaD',0.1,'alphaF',0); %alpha-values for detection statistical tests
        detectionParam.visual = 1; %1 to see image with detected features, 0 otherwise
        detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
        detectionParam.bitDepth = 16; %Camera bit depth
        detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
        detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
        detectionParam.integWindow = 0; %number of frames before and after a frame for time integration
        
%         %background info ...
%         background.imageDir = bgImageDir{iMovie};
%         background.filenameBase = bgFilenameBase{iMovie};
%         background.alphaLocMaxAbs = 1e-15;
%         detectionParam.background = background;
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionFyn.mat'; %name of file where input and output are saved
        % saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
