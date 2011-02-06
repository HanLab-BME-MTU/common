classdef TrackingProcess < Process
% A class definition for a generic tracking process.
%
% Chuangang Ren, 11/2010

properties(SetAccess = protected, GetAccess = public)
    
    
    channelIndex_ % The index of channel to process
    filename_ % file name of result data
    overwrite_ = 0; % If overwrite the original result MAT file
    
end

methods(Access = public)
    
    function obj = TrackingProcess(owner, outputDir, channelIndex, funParams )
       
        if nargin == 0
            super_args = {};
        else
            super_args{1} = owner;
            super_args{2} = 'Sub-Resolution';
        end
        obj = obj@Process(super_args{:});
        
        if nargin < 2 || isempty(outputDir)
            outputDir = owner.outputDirectory_ ;
        end   
        
        if nargin < 3 || isempty(channelIndex) % Default channel Index
            channelIndex = 1:length(owner.channels_);
        end        
        
        if nargin < 4 || isempty(funParams)
            
            % --------------- gapCloseParam ----------------
            
            funParams.gapCloseParam.timeWindow = 5; %IMPORTANT maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; % (SORT OF FLAG: 4 options for user) 1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.
            funParams.gapCloseParam.diagnostics = 1; %FLAG 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.            
            
            % --------------- costMatrices ----------------
            
            % Linking:
            funParams.costMatrices(1).funcName = 'costMatLinearMotionLink2';

            parameters.linearMotion = 0; %FLAG use linear motion Kalman filter.
            parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            parameters.maxSearchRadius = 5; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            parameters.nnWindow = funParams.gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
            parameters.diagnostics = [2 obj.owner_.nFrames_-1]; %FLAG THEN NUMBERS if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

            funParams.costMatrices(1).parameters = parameters;
            clear parameters
            
            % Gap Closing:
            funParams.costMatrices(2).funcName = 'costMatLinearMotionCloseGaps2';
            parameters.linearMotion = 0; %use linear motion Kalman filter.

            parameters.minSearchRadius = 2; %minimum allowed search radius.
            parameters.maxSearchRadius = 5; %maximum allowed search radius.
            parameters.brownStdMult = 3*ones(funParams.gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

            parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            parameters.nnWindow = funParams.gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
            parameters.brownScaling = [0.5 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
            parameters.timeReachConfB = 5; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            parameters.ampRatioLimit = [0.5 2]; % (FLAG + VALUES small-big value) for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

            % If parameters.linearMotion = 1
            parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
            parameters.linStdMult = 3*ones(funParams.gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
            parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            parameters.timeReachConfL = 5; %similar to timeReachConfB, but for the linear part of the motion.
            parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
            % ---------------------------------
            
            parameters.gapPenalty = 1; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
            parameters.resLimit = []; % text field resolution limit, which is generally equal to 3 * point spread function sigma.

            funParams.costMatrices(2).parameters = parameters;
            clear parameters            
            
            % --------------- kalmanFunctions ----------------
            
            funParams.kalmanFunctions.reserveMem  = 'kalmanResMemLM';
            funParams.kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
            funParams.kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
            funParams.kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
            
            % --------------- saveResults ----------------
            
            funParams.saveResults.dir = [outputDir  filesep 'Tracking' filesep]; %directory where to save input and output
            funParams.saveResults.filename = []; % Note: channel-specific
            
            % --------------- Others ----------------
            
            funParams.verbose = 1;
            funParams.probDim = 2;
            
        end
        
        obj.funParams_ = funParams;
        
        obj.outParams_ = cell(1, length(owner.channels_));
        obj.channelIndex_ = channelIndex;
        obj.filename_ = 'tracking_result.mat';
        
        obj.funName_ = @trackCloseGapsKalmanSparse;
        
        % ---------------- Visualization Parameters --------------------
        
        % Tool 1: plotTracks2D
        
        obj.visualParams_.pt2D.timeRange = [1 owner.nFrames_];
        obj.visualParams_.pt2D.colorTime = 'k';
        obj.visualParams_.pt2D.markerType = 'none';
        obj.visualParams_.pt2D.indicateSE = 1;
        obj.visualParams_.pt2D.newFigure = 1;
        obj.visualParams_.pt2D.image = [];
        obj.visualParams_.pt2D.imageDir = []; % Not in original function
        obj.visualParams_.pt2D.flipXY = 0;
        obj.visualParams_.pt2D.ask4sel = 1;
        obj.visualParams_.pt2D.offset = [0 0];
        
        % Tool 2: plotCompTrack
        
        obj.visualParams_.pct.trackid = 1; % Not in original function
        obj.visualParams_.pct.plotX = 1;
        obj.visualParams_.pct.plotY = 1;
        obj.visualParams_.pct.plotA = 1;
        obj.visualParams_.pct.inOneFigure = 1;
        obj.visualParams_.pct.plotAggregState = 1;
        
        % Tool 3: overlayTracksMovieNew
        
        obj.visualParams_.otmn.startend = [1 owner.nFrames_];
        obj.visualParams_.otmn.dragtailLength = 5;
        obj.visualParams_.otmn.saveMovie = 0;
        obj.visualParams_.otmn.movieName = [];
        obj.visualParams_.otmn.dir2saveMovie = funParams.saveResults.dir;
        obj.visualParams_.otmn.filterSigma = 0;
        obj.visualParams_.otmn.classifyGaps = 1;
        obj.visualParams_.otmn.highlightES = 1;
        obj.visualParams_.otmn.showRaw = 0;
        obj.visualParams_.otmn.imageRange = []; % TO DO in GUI
        obj.visualParams_.otmn.onlyTracks = 0;
        obj.visualParams_.otmn.classifyLft = 1;
        obj.visualParams_.otmn.diffAnalysisRes = [];
        obj.visualParams_.otmn.intensityScale = 1;
        obj.visualParams_.otmn.colorTracks = 0;
        file = owner.getImageFileNames(1);
        obj.visualParams_.otmn.firstImageFile = [owner.channels_(1).channelPath_ filesep file{1}{1}];        
    end
    
    function sanityCheck(obj) 
    end
    

    
    function setChannelIndex(obj, index)
        % Set channel index
        if any(index > length(obj.owner_.channels_))
           error ('User-defined: channel index is larger than the number of channels.') 
        end
        obj.channelIndex_ = index;
    end    
    
    % Set result file name
    function setFileName(obj, name)
       obj.filename_ = name;
    end   
    
    % Set overwrite
    function setOverwrite (obj, i)
       obj.overwrite_ = i;
    end
    
    function runProcess(obj)
    % Run the process!
        iDetection = find(cellfun(@(x)(isa(x,'DetectionProcess')), obj.owner_.processes_),1);
        
        if isempty(iDetection)
           error('The detection step has not been set up yet.') 
        end
        
        if any(cellfun(@(x)isempty(x), obj.owner_.processes_{iDetection}.outParams_(obj.channelIndex_)))
           error('One or more channels specified has not been processed by detection step.') 
        end
        
        
        for i = obj.channelIndex_
            movieInfo = obj.owner_.processes_{iDetection}.outParams_{i}.movieInfo;
            obj.funParams_.saveResults.filename = ['Channel_' num2str(i) '_' obj.filename_];
            
            %Check/create directory
            if ~exist(obj.funParams_.saveResults.dir,'dir')
                mkdir(obj.funParams_.saveResults.dir)
            end
            
            if ~obj.overwrite_
                % file name enumaration
                obj.funParams_.saveResults.filename = enumFileName(obj.funParams_.saveResults.dir, obj.funParams_.saveResults.filename);
            end
            
            % Test (commentable)
%             movieInfo, obj.funParams_.costMatrices, obj.funParams_.gapCloseParam, ...
%                 obj.funParams_.kalmanFunctions, obj.funParams_.probDim, obj.funParams_.saveResults, obj.funParams_.verbose
            
            % Call function
            [obj.outParams_{i}.tracksFinal, obj.outParams_{i}.kalmanInfoLink, obj.outParams_{i}.errFlag] = ...
                obj.funName_(movieInfo, obj.funParams_.costMatrices, obj.funParams_.gapCloseParam, ...
                obj.funParams_.kalmanFunctions, obj.funParams_.probDim, obj.funParams_.saveResults, obj.funParams_.verbose);
            
        end
        
%         [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
%             costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
        
    end

    
end

methods (Static)
    
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end
            description = 'This is help of tracking.';            
            paramList = {''};
                         
            paramDesc = {''};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end   
    
end


end