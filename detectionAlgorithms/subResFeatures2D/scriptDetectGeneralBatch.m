
%% define batch job locations

%image locations
imageDir = {...
    '/home/kj35/tmpData/CD36/2010-12-10-CD36andMTs/2010-12-09_CD36Tub_export/Coverslip02/movie47/images/',...
    };

%file name bases
filenameBase = {...
    'movie_47_',...
    };

%directory for saving results
saveResDir = {...
    '/home/kj35/tmpData/CD36/2010-12-10-CD36andMTs/2010-12-09_CD36Tub_export/Coverslip02/movie47/analysis/',...
    };

%background image locations
% bgImageDir = {...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_1/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_5/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_6/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_10/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_12/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_17/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_18/bgImages/',...
%     '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_20/bgImages/',...
%     };

%background file name bases
% bgFilenameBase = {...
%     'crop_WT1_1_D3D_',...
%     'crop_WT1_5_D3D_',...
%     'crop_WT1_6_D3D_',...
%     'crop_WT1_10_D3D_',...
%     'crop_WT1_12_D3D_',...
%     'crop_WT1_17_D3D_',...
%     'crop_WT1_18_D3D_',...
%     'crop_WT1_20_D3D_',...
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
        movieParam.lastImageNum = 443; %number of last image in movie
        movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).
        
        %% detection parameters
        detectionParam.psfSigma = 1.5; %point spread function sigma (in pixels)
        detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
        detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
        detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
        detectionParam.bitDepth = 16; %Camera bit depth
        detectionParam.alphaLocMax = 0.01; %alpha-value for initial detection of local maxima
        detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
        detectionParam.integWindow = 1; %number of frames before and after a frame for time integration
        
        %background info ...
        %         background.imageDir = bgImageDir{iMovie};
        %         background.filenameBase = bgFilenameBase{iMovie};
        %         background.alphaLocMaxAbs = 1e-15;
        %         detectionParam.background = background;
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detectionAll1.mat'; %name of file where input and output are saved
        % saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
