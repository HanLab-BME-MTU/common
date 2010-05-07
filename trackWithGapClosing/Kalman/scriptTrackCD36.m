
%% START PARAMETERS TO VARY

%where to save results
saveResults.dir = '???\CD36\con_cy3_02\analysis\'; %directory where to save input and output
saveResults.filename = 'tracksAll_1.mat'; %name of file where input and output are saved

%flag for linear motion modeling
%1 = use linear motion Kalman filter. 0 = don't use.
parameters.linearMotion = 0;

%maximum allowed search radius
parameters.maxSearchRadius = 2; 

%diagnostics: to plot the histogram of linking distances up to a certain
%frame, indicate its number
parameters.diagnostics = 98;

%gap closing time window
gapCloseParam.timeWindow = 7;

%diagnostics: 1 to plot a histogram of gap lengths in the end; 0 or empty
%otherwise
gapCloseParam.diagnostics = 0; 

%% END PARAMETERS TO VARY






%% other general gap closing parameters

gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.

%% other parameters for cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatLinearMotionLink2';

%parameters

parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
% parameters.nnWindow = 10; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam = []; %Kalman filter initialization parameters.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps2';

%parameters

parameters.linearMotion = costMatrices(1).parameters.linearMotion; %1 = use linear motion Kalman filter. 0 = don't use.

parameters.maxSearchRadius = costMatrices(1).parameters.maxSearchRadius; %maximum allowed search radius.

parameters.minSearchRadius = 2; %minimum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
parameters.timeReachConfB = gapCloseParam.timeWindow; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).

parameters.ampRatioLimit = [0 Inf]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
% parameters.nnWindow = 10; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
parameters.timeReachConfL = gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.1; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% tracking function call

[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

