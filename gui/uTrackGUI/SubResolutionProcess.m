classdef SubResolutionProcess < DetectionProcess
    % A subclass of the detection process
    % Chuangang Ren
    % 11/2010
    
    properties(SetAccess = protected, GetAccess = public)
        
        filenameBase_   % cell array storing the base of file name
        digits4Enum_    % cell array stroring the number of digits for frame enumeration (1-4)
        filename_       % file name of result data
    end
    
    methods (Access = public)
        function obj = SubResolutionProcess(owner, outputDir, channelIndex, funParams )
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = SubResolutionProcess.getName;
            super_args{3} = @detectSubResFeatures2D_StandAlone;
            
            if nargin < 2 || isempty(outputDir)
                outputDir = owner.outputDirectory_ ; % package folder
            end
            
            if nargin < 3 || isempty(channelIndex) % Default channel Index
                channelIndex = 1:length(owner.channels_);
            end
            
            super_args{4} = channelIndex;
            
            
            if nargin < 4 || isempty(funParams)  % Default funParams
                funParams = SubResolutionProcess.getDefaultParams(owner,outputDir);
            end
            
            super_args{5} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
            
            % Visual parameters ( Default: channel 1 )
            obj.visualParams_.startend = [1 owner.nFrames_];
            obj.visualParams_.saveMovie = 1;
            obj.visualParams_.movieName = [];
            obj.visualParams_.dir2saveMovie = funParams.saveResults.dir;
            obj.visualParams_.filterSigma = 0;
            obj.visualParams_.showRaw = 1;
            obj.visualParams_.intensityScale = 1;
            file = owner.getImageFileNames(1);
            obj.visualParams_.firstImageFile = [owner.channels_(1).channelPath_ filesep file{1}{1}];
            
            
            % Get file name base and digits for enumeration
            [obj.filenameBase_ obj.digits4Enum_] = SubResolutionProcess.getFilenameBody(owner);
            obj.filename_ = 'detection_result.mat';
        end
        
        
        function setFileName(obj, name)
            % Set result file name
            obj.filename_ = name;
        end
        
        function OK = checkChannelOutput(obj,iChan)
            
            %Checks if the selected channels have valid output files
            nChanTot = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan)
                iChan = 1:nChanTot;
            end
            %Makes sure there's at least one .mat file in the speified
            %directory
            OK =  arrayfun(@(x)(x <= nChanTot && ...
                x > 0 && isequal(round(x),x) && ...
                exist(obj.outFilePaths_{x},'file')),iChan);
        end
        
        function run(obj)
            % Run the process!
            obj.success_=false;
            
            for i = obj.channelIndex_
                
                obj.funParams_.movieParam.imageDir = [obj.owner_.channels_(i).channelPath_ filesep];
                obj.funParams_.movieParam.filenameBase = obj.filenameBase_{i};
                obj.funParams_.movieParam.digits4Enum = obj.digits4Enum_{i};
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
                %                 obj.funParams_.movieParam,obj.funParams_.detectionParam,obj.funParams_.saveResults
                obj.funName_(obj.funParams_.movieParam, obj.funParams_.detectionParam, obj.funParams_.saveResults);
                obj.setOutFilePath(i,[obj.funParams_.saveResults.dir filesep obj.funParams_.saveResults.filename]);
            end
            obj.success_=true;
            obj.procChanged_=false;
            obj.setDateTime;
            obj.owner_.save;
        end
        function hfigure = resultDisplay(obj,fig,procID)
            % Display the output of the process
            
            % Copied and pasted from the old uTrackPackageGUI
            % but there is definitely some optimization to do
            % Check for movie output before loading the GUI
            chan = [];
            for i = 1:length(obj.owner_.channels_)
                if obj.checkChannelOutput(i)
                    chan = i;
                    break
                end
            end
            
            if isempty(chan)
                warndlg('The current step does not have any output yet.','No Output','modal');
                return
            end
            
            % Make sure detection output is valid
            load(obj.outFilePaths_{chan},'movieInfo');
            firstframe = [];
            for i = 1:length(movieInfo)
                
                if ~isempty(movieInfo(i).amp)
                    firstframe = i;
                    break
                end
            end
            
            if isempty(firstframe)
                warndlg('The detection result is empty. There is nothing to visualize.','Empty Output','modal');
                return
            end
            
            if isa(obj, 'Process')
                hfigure = detectionVisualGUI('mainFig', fig, procID);
            else
                error('User-defined: the input is not a Process object.')
            end
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'movieInfo'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'SubResolutionProcess'));
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{iChan},output{:});
           
            if numel(ip.Results.iFrame)>1,
                varargout{1}=s.(output{1});
            else
                varargout{1}=s.(output{1})(iFrame);
            end
        end
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Sub-resolution objects';
            output(1).var='movieInfo';
            output(1).formatData=@(x) horzcat(x.yCoord(:,1),x.xCoord(:,1));
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none','Color',colors(x,:));
        end
        
    end
    methods (Static)
        
        function [base digits4Enum]= getFilenameBody(owner)
            % Get the base of file name for all channels in movie data
            % "owner"
            
            fileNames = owner.getImageFileNames;
            base = cell(1, length(owner.channels_));
            digits4Enum = cell(1, length(owner.channels_));
            
            for i = 1 : length(owner.channels_)
                
                [x1 base{i} digits4Enum{i} x4] = getFilenameBody(fileNames{i}{1});
                digits4Enum{i} = length(digits4Enum{i});
            end
            
        end
        function name = getName()
            name = 'Sub-Resolution Detection';
        end
        function h = GUI()
            h = subResolutionProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % movieParam
            funParams.movieParam.imageDir = owner.channels_(1).channelPath_; % Note: channel-specific
            funParams.movieParam.filenameBase = []; % Note: channel-specific
            funParams.movieParam.firstImageNum = 1;
            funParams.movieParam.lastImageNum = owner.nFrames_;
            funParams.movieParam.digits4Enum = []; % Note: channel-specific
            
            % detectionParam
            %                 funParams.detectionParam.psfSigma = [];
            %                 funParams.detectionParam.bitDepth = owner.camBitdepth_;
            funParams.detectionParam.alphaLocMax = .05;
            funParams.detectionParam.integWindow = 0;
            funParams.detectionParam.doMMF = 0;
            funParams.detectionParam.testAlpha = struct('alphaR', .05,'alphaA', .05, 'alphaD', .05,'alphaF',0);
            funParams.detectionParam.numSigmaIter = 0;
            funParams.detectionParam.visual = 0;
            funParams.detectionParam.background = [];
            
            % saveResults
            %                 funParams.OutputDirectory = [outputDir  filesep 'Sub_Resolution_Detection' filesep];
            funParams.saveResults.dir = [outputDir  filesep 'Sub_Resolution_Detection' filesep];
            funParams.saveResults.filename = []; % Note: channel-specific
            
            % Set up psfSigma and bitDepth
            na = owner.numAperture_;
            ps = owner.pixelSize_;
            wl = owner.channels_(1).emissionWavelength_;
            bd = owner.camBitdepth_;
            
            if ~isempty( na ) && ~isempty( ps ) && ~isempty( wl )
                funParams.detectionParam.psfSigma = 0.21*wl/na/ps;
            else
                funParams.detectionParam.psfSigma = [];
            end
            
            if ~isempty(bd)
                funParams.detectionParam.bitDepth = bd;
            else
                funParams.detectionParam.bitDepth = [];
            end
            
        end

    end
    
end