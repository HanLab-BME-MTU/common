
function batchJob = setupBatchDetectTrack(numJobs)

for iJob = 1 : numJobs
      
    disp(['---Please enter in the file the information and parameters for movie # ' num2str(iJob) ','])
    input('    then save the file and press enter to proceed to the next movie');
    
    %% things that have to be changed ...
    
    %movie information
    movieParam.imageDir = 'Z:\Shared\Waterman\IntegrinGroup\090707\integrin fn2 mn2 hr36\cell A\cell A 50ms 150x B\'; %directory where images are
    movieParam.filenameBase = 'cell A 50ms 150x B'; %image file name base
    movieParam.firstImageNum = 1; %number of first image in movie
    movieParam.lastImageNum = 10; %number of last image in movie
    movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).
    
    %save detection results
    saveResultsDetect.dir = 'Z:\Shared\Waterman\IntegrinGroup\090707\integrin fn2 mn2 hr36\cell A\cell A 50ms 150x B'; %directory where to save input and output
    saveResultsDetect.filename = 'testDetection1.mat'; %name of file where input and output are saved
    
    %save tracking results
    saveResultsTrack.dir = 'Z:\Shared\Waterman\IntegrinGroup\090707\integrin fn2 mn2 hr36\cell A\cell A 50ms 150x B'; %directory where to save input and output
    saveResultsTrack.filename = 'testTracking1.mat'; %name of file where input and output are saved

    %% things that can change ...
    
    %detection parameters
    detectionParam.psfSigma = 0.8; %point spread function sigma (in pixels)
    detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
    detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
    detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
    detectionParam.bitDepth = 16; %Camera bit depth
    detectionParam.alphaLocMax = 0.02; %alpha-value for initial detection of local maxima
    detectionParam.numSigmaIter = 10; %maximum number of iterations for PSF sigma estimation
    detectionParam.integWindow = 1; %number of frames before and after a frame for time integration
    
    %general gap closing parameters
    gapCloseParam.timeWindow = 10; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
    gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
    gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
    
    %optional input:
    gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
    
    %cost matrix for frame-to-frame linking
    
    %function name
    costMatrices(1).funcName = 'costMatLinearMotionLink';
    
    %parameters
    
    parameters.linearMotion = 1; %use linear motion Kalman filter.
    
    parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
    parameters.maxSearchRadius = 5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
    parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
    
    parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
    parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
    
    parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
    
    %optional input
    parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.
    
    costMatrices(1).parameters = parameters;
    clear parameters
    
    %cost matrix for gap closing
    
    %function name
    costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';
    
    %parameters
    
    %needed all the time
    parameters.linearMotion = 1; %use linear motion Kalman filter.
    
    parameters.minSearchRadius = 2; %minimum allowed search radius.
    parameters.maxSearchRadius = 5; %maximum allowed search radius.
    parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
    parameters.timeReachConfB = 3; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
    
    parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
    
    parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
    
    parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
    parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
    
    parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
    parameters.timeReachConfL = gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
    parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
    
    %optional; if not input, 1 will be used (i.e. no penalty)
    parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
    
    %optional; to calculate MS search radius
    %if not input, MS search radius will be the same as gap closing search radius
    parameters.resLimit = 5.1; %resolution limit, which is generally equal to 3 * point spread function sigma.
    
    costMatrices(2).parameters = parameters;
    clear parameters
    
    %Kalman filter function names
    
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
    
    %verbose
    verbose = 1;
    
    %problem dimension
    probDim = 2;
    
    batchJob(iJob).movieParam = movieParam;
    batchJob(iJob).detectionParam = detectionParam;
    batchJob(iJob).saveResultsDetect = saveResultsDetect;
    batchJob(iJob).gapCloseParam = gapCloseParam;
    batchJob(iJob).costMatrices = costMatrices;
    batchJob(iJob).kalmanFunctions = kalmanFunctions;
    batchJob(iJob).verbose = verbose;
    batchJob(iJob).probDim = probDim;
    batchJob(iJob).saveResultsTrack = saveResultsTrack;
    
end
