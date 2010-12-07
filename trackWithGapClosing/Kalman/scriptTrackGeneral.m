
%% general gap closing parameters
gapCloseParam.timeWindow = 6; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatLinearMotionLink2';

%parameters

parameters.linearMotion = 1; %use linear motion Kalman filter.

parameters.minSearchRadius = 3; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = 3; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

% parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
parameters.kalmanInitParam.initVelocity = orientVec; %Kalman filter initialization parameters.

%optional input
parameters.diagnostics = [199]; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps2';

%parameters

%needed all the time
parameters.linearMotion = 1; %use linear motion Kalman filter.

parameters.minSearchRadius = 3; %minimum allowed search radius.
parameters.maxSearchRadius = 3; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

parameters.brownScaling = [0.5 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
parameters.timeReachConfB = 1; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

parameters.ampRatioLimit = [0 Inf]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.25; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotionDiscont';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%saveResults
saveResults.dir = '/home/kj35/tmpData/Martin/2010_10_19_WT/WT1_1/analysis/'; %directory where to save input and output
saveResults.filename = 'trackingIllustrate.mat'; %name of file where input and output are saved
% saveResults = 0; %don't save results

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% tracking function call

[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

% saveResults.filename = 'tracks1Detection2_Frames0001to1200.mat';
% [tracksFinal01,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(1:1200),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2_Frames1201to2400.mat';
% [tracksFinal02,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(1201:2400),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2_Frames2401to3600.mat';
% [tracksFinal03,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(2401:3600),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1DetectionMF1Frames3601to4800.mat';
% [tracksFinal04,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(3601:4800),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1DetectionMF1Frames4801to6000.mat';
% [tracksFinal05,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(4801:6000),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames6001to7200.mat';
% [tracksFinal06,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(6001:7200),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames7201to8200.mat';
% [tracksFinal07,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(7201:8200),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames4201to4800.mat';
% [tracksFinal08,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(4201:4800),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames4801to5400.mat';
% [tracksFinal09,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(4801:5400),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames5401to6000.mat';
% [tracksFinal10,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(5401:6000),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames6001to6600.mat';
% [tracksFinal11,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(6001:6600),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% 
% saveResults.filename = 'tracks1Detection2Frames6601to7200.mat';
% [tracksFinal12,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(6601:7200),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

