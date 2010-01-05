
%% define batch job locations

%image locations
imageDir = {'/orchestra/groups/lccb-receptors/Hiro/091208Cy3CD36jasp/expo2/con14/',...
    '/orchestra/groups/lccb-receptors/Hiro/091208Cy3CD36jasp/expo2/con16/'};

%file name bases
filenameBase = {'0204Fab36_SPT_jasp_',...
    '0204Fab36_SPT_jasp_'};

%directory for saving results
saveResDir = {'/orchestra/groups/lccb-receptors/Hiro/091208Cy3CD36jasp/expo2/con14/',...
    '/orchestra/groups/lccb-receptors/Hiro/091208Cy3CD36jasp/expo2/con16/'};

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
        movieParam.lastImageNum = 50; %number of last image in movie
        movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).
        
        %% detection parameters
        detectionParam.psfSigma = 1.7; %point spread function sigma (in pixels)
        detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
        detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
        detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
        detectionParam.bitDepth = 16; %Camera bit depth
        detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
        detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
        detectionParam.integWindow = 1; %number of frames before and after a frame for time integration
        
        %% save results
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'detection1.mat'; %name of file where input and output are saved
        % saveResults = 0;
        
        %% run the detection function
        [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
            detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
