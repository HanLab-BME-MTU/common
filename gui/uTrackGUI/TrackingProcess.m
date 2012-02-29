classdef TrackingProcess < DataProcessingProcess
    % A class definition for a generic tracking process.
    %
    % Chuangang Ren, 11/2010
    % Sebastien Besson (last modified Dec 2011)

   
    methods(Access = public)
        
        function obj = TrackingProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = TrackingProcess.getName;
                super_args{3} = @trackMovie;
                if isempty(funParams)
                    funParams = TrackingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
            
            % ---------------- Visualization Parameters --------------------
            
            % Tool 1: plotTracks2D
            
            obj.visualParams_.pt2D.timeRange = [1 owner.nFrames_];
            obj.visualParams_.pt2D.colorTime = '3';
            obj.visualParams_.pt2D.markerType = 'none';
            obj.visualParams_.pt2D.indicateSE = 0;
            obj.visualParams_.pt2D.newFigure = 1;
            obj.visualParams_.pt2D.image = [];
            obj.visualParams_.pt2D.imageDir = []; % Not in original function
            obj.visualParams_.pt2D.flipXY = 0;
            obj.visualParams_.pt2D.ask4sel = 0;
            obj.visualParams_.pt2D.offset = [0 0];
            obj.visualParams_.pt2D.minLength = 1;
            
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
            obj.visualParams_.otmn.saveMovie = 1;
            obj.visualParams_.otmn.movieName = [];
            obj.visualParams_.otmn.dir2saveMovie = funParams.OutputDirectory;
            obj.visualParams_.otmn.filterSigma = 0;
            obj.visualParams_.otmn.classifyGaps = 0;
            obj.visualParams_.otmn.highlightES = 0;
            obj.visualParams_.otmn.showRaw = 1;
            obj.visualParams_.otmn.imageRange = []; % TO DO in GUI
            obj.visualParams_.otmn.onlyTracks = 0;
            obj.visualParams_.otmn.classifyLft = 0;
            obj.visualParams_.otmn.diffAnalysisRes = [];
            obj.visualParams_.otmn.intensityScale = 1;
            obj.visualParams_.otmn.colorTracks = 1;
            obj.visualParams_.otmn.minLength = 1;
            file = owner.getImageFileNames(1);
            obj.visualParams_.otmn.firstImageFile = [owner.channels_(1).channelPath_ filesep file{1}{1}];
            
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'tracksFinal'};
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',[],@(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{iChan},output{:});
            tracksFinal=s.tracksFinal;
            
            if ~isempty(iFrame),
                startTime = arrayfun(@(x) x.seqOfEvents(1,1),tracksFinal);
                endTime = arrayfun(@(x) x.seqOfEvents(2,1),tracksFinal);
                isValid = (iFrame>=startTime &iFrame<=endTime);
                data.x = arrayfun(@(x,y) x.tracksCoordAmpCG(1,1:8:1+8*(iFrame-y)),...
                    tracksFinal(isValid),startTime(isValid),'Unif',0);
                data.y = arrayfun(@(x,y) x.tracksCoordAmpCG(1,2:8:2+8*(iFrame-y)),...
                    tracksFinal(isValid),startTime(isValid),'Unif',0);
                data.label = find(isValid);
                varargout{1}=data;
            else
                varargout{1} = tracksFinal;
            end
            
        end
        
        
        
        function hfigure = resultDisplay(obj,varargin)
            % Display the output of the process
            if ~isempty(obj.getPackage) && ...
                    isa(obj.owner_.packages_{obj.getPackage},'UTrackPackage')
                % Check for movie output before loading the GUI
                iChan = find(obj.checkChannelOutput,1);
                if isempty(iChan)
                    warndlg('The current step does not have any output yet.','No Output','modal');
                    return
                end
                
                % Make sure detection output is valid
                tracksFinal=obj.loadChannelOutput(iChan);
                if isempty(tracksFinal)
                    warndlg('The tracking result is empty. There is nothing to visualize.','Empty Output','modal');
                    return
                end
                
                hfigure = trackingVisualGUI('mainFig', varargin{:});
            else
                hfigure=resultDisplay@Process(obj,varargin{:});
            end
        end
        
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Tracks';
            output(1).var='tracksFinal';
            output(1).formatData=[];
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)TracksDisplay('Color',colors(x,:));
        end
        
        
    end
    methods(Static)
        function name = getName()
            name = 'Tracking';
        end
        function h = GUI()
            h= @trackingProcessGUI;
        end

        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
                       
            funParams.ChannelIndex =1:numel(owner.channels_);
            funParams.DetProcessIndex = [];
            funParams.OutputDirectory = [outputDir  filesep 'Tracking'];
            
            % --------------- gapCloseParam ----------------
            funParams.gapCloseParam.timeWindow = 5; %IMPORTANT maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; % (SORT OF FLAG: 4 options for user) 1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.
            funParams.gapCloseParam.diagnostics = 1; %FLAG 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            
            
            % --------------- kalmanFunctions ----------------
            
            funParams.kalmanFunctions.reserveMem  = TrackingProcess.getKalmanReserveMemFunctions(1).funcName;
            funParams.kalmanFunctions.initialize  = TrackingProcess.getKalmanInitializeFunctions(1).funcName;
            funParams.kalmanFunctions.calcGain    = TrackingProcess.getKalmanCalcGainFunctions(1).funcName;
            funParams.kalmanFunctions.timeReverse = TrackingProcess.getKalmanTimeReverseFunctions(1).funcName;
            
            
            % --------------- saveResults ----------------
            funParams.saveResults.export = 0; %FLAG allow additional export of the tracking results into matrix
            
            % --------------- Others ----------------
            
            funParams.verbose = 1;
            funParams.probDim = 2;
            
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            

        end
        
        function functions = getKalmanReserveMemFunctions(varargin)
            functions(1).name = 'Brownian + Directed motion models';
            functions(1).funcName = func2str(@kalmanResMemLM);
            
            ip=inputParser;
            ip.addOptional('index',1:length(functions),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            functions=functions(index);      
        end
        
        function functions = getKalmanInitializeFunctions(varargin)
            functions(1).name = 'Brownian + Directed motion models';
            functions(1).funcName = func2str(@kalmanInitLinearMotion);
            functions(1).GUI=@kalmanInitializationGUI;
            functions(2).name = 'Microtubule dynamics';
            functions(2).funcName = func2str(@plusTipKalmanInitLinearMotion);
            functions(2).GUI=@kalmanInitializationGUI;
            ip=inputParser;
            ip.addOptional('index',1:length(functions),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            functions=functions(index);      
        end
        
        function functions = getKalmanCalcGainFunctions(varargin)
            functions(1).name = 'Brownian + Directed motion models';
            functions(1).funcName = func2str(@kalmanGainLinearMotion);
            functions(2).name = 'Microtubule dynamics';
            functions(2).funcName = func2str(@plusTipKalmanGainLinearMotion);
            ip=inputParser;
            ip.addOptional('index',1:length(functions),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            functions=functions(index);      
        end
        
        function functions = getKalmanTimeReverseFunctions(varargin)
            functions(1).name = 'Brownian + Directed motion models';
            functions(1).funcName = func2str(@kalmanReverseLinearMotion);
            
            ip=inputParser;
            ip.addOptional('index',1:length(functions),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            functions=functions(index);      
        end     
        
        
        function costMatrix = getDefaultLinkingCostMatrices(owner,timeWindow,varargin)
            
            % Brownian + Directed motion models
            costMatrices(1).name = 'Brownian + Directed motion models';
            costMatrices(1).funcName = func2str(@costMatRandomDirectedSwitchingMotionLink);
            costMatrices(1).GUI = @costMatRandomDirectedSwitchingMotionLinkGUI;
            costMatrices(1).parameters.linearMotion = 0; % use linear motion Kalman filter.
            costMatrices(1).parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            costMatrices(1).parameters.maxSearchRadius = 5; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            costMatrices(1).parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            costMatrices(1).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            costMatrices(1).parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            costMatrices(1).parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
            costMatrices(1).parameters.diagnostics = owner.nFrames_-1;
            
            % plusTip markers
            plusTipCostMatrix.name = 'Microtubule dynamics';
            plusTipCostMatrix.funcName = func2str(@plusTipCostMatLinearMotionLink);
            plusTipCostMatrix.GUI = @plusTipCostMatLinearMotionLinkGUI;
            plusTipCostMatrix.parameters.linearMotion = 1; % use linear motion Kalman filter.
            plusTipCostMatrix.parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            plusTipCostMatrix.parameters.maxSearchRadius = 10; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            plusTipCostMatrix.parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            plusTipCostMatrix.parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            plusTipCostMatrix.parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            plusTipCostMatrix.parameters.kalmanInitParam.initVelocity = []; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.kalmanInitParam.convergePoint = []; %Kalman filter initialization parameters.  
            plusTipCostMatrix.parameters.kalmanInitParam.searchRadiusFirstIteration = 20; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.diagnostics = [];
            costMatrices(2)=plusTipCostMatrix;
            
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addRequired('timeWindow',@isscalar);
            ip.addOptional('index',1:length(costMatrices),@isvector);
            ip.parse(owner,timeWindow,varargin{:});
            index = ip.Results.index;
            costMatrix=costMatrices(index);          
        end
        
        function costMatrix = getDefaultGapClosingCostMatrices(owner,timeWindow,varargin)
            
            % Linear motion
            costMatrices(1).name = 'Brownian + Directed motion models';
            costMatrices(1).funcName = func2str(@costMatRandomDirectedSwitchingMotionCloseGaps);
            costMatrices(1).GUI = @costMatRandomDirectedSwitchingMotionCloseGapsGUI;
            costMatrices(1).parameters.linearMotion = 0; %use linear motion Kalman filter.
            
            costMatrices(1).parameters.minSearchRadius = 2; %minimum allowed search radius.
            costMatrices(1).parameters.maxSearchRadius = 5; %maximum allowed search radius.
            costMatrices(1).parameters.brownStdMult = 3*ones(timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
            
            costMatrices(1).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            costMatrices(1).parameters.nnWindow = timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
            costMatrices(1).parameters.brownScaling = [0.5 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
            costMatrices(1).parameters.timeReachConfB = timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            costMatrices(1).parameters.ampRatioLimit = [0.5 2]; % (FLAG + VALUES small-big value) for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
            
            % If parameters.linearMotion = 1
            costMatrices(1).parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
            costMatrices(1).parameters.linStdMult = 3*ones(timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
            costMatrices(1).parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            costMatrices(1).parameters.timeReachConfL = timeWindow; %similar to timeReachConfB, but for the linear part of the motion.
            costMatrices(1).parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
            % ---------------------------------
            
            costMatrices(1).parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
            costMatrices(1).parameters.resLimit = []; % text field resolution limit, which is generally equal to 3 * point spread function sigma.
            
            % Linear motion
            plusTipCostMatrix.name = 'Microtubule dynamics';
            plusTipCostMatrix.funcName = func2str(@plusTipCostMatCloseGaps);
            plusTipCostMatrix.GUI = @plusTipCostMatCloseGapsGUI;
            plusTipCostMatrix.parameters.maxFAngle = 30; %use linear motion Kalman filter.
            plusTipCostMatrix.parameters.maxBAngle = 10; %use linear motion Kalman filter.
            plusTipCostMatrix.parameters.backVelMultFactor = 1.5;
            plusTipCostMatrix.parameters.fluctRad = 1.0;
            plusTipCostMatrix.parameters.breakNonLinearTracks = false;
            costMatrices(2)=plusTipCostMatrix;
            
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addRequired('timeWindow',@isscalar);
            ip.addOptional('index',1:length(costMatrices),@isvector);
            ip.parse(owner,timeWindow,varargin{:});
            index = ip.Results.index;
            costMatrix=costMatrices(index);      
        end    
        
    end
end