classdef SubResolutionProcess < DetectionProcess
    % A concrete class for detecting objects using Gaussian mixture-model fitting
    % Chuangang Ren 11/2010
    % Sebastien Besson (last modified Dec 2011)
    
    methods (Access = public)
        function obj = SubResolutionProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            outputDir = ip.Results.outputDir;
            funParams = ip.Results.funParams;
            
            
            % Constructor of the SubResolutionProcess
            
            super_args{1} = owner;
            super_args{2} = SubResolutionProcess.getName;
            super_args{3} = @detectMovieSubResFeatures;
            if isempty(funParams)  % Default funParams
                funParams = SubResolutionProcess.getDefaultParams(owner,outputDir);
            end
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
            
            % Visual parameters ( Default: channel 1 )
            obj.visualParams_.startend = [1 owner.nFrames_];
            obj.visualParams_.saveMovie = 1;
            obj.visualParams_.movieName = [];
            obj.visualParams_.dir2saveMovie = funParams.OutputDirectory;
            obj.visualParams_.filterSigma = 0;
            obj.visualParams_.showRaw = 1;
            obj.visualParams_.intensityScale = 1;
            if isa(owner.channels_(1).channelPath_, 'omero.model.Pixels')
                obj.visualParams_.firstImageFile = owner.getChannelPaths{1};
            elseif exist(owner.getChannelPaths{1}, 'file')==2  
                obj.visualParams_.firstImageFile = owner.getChannelPaths{1};
            else
                obj.visualParams_.firstImageFile = [owner.getChannelPaths{1} filesep owner.getImageFileNames{1}{1}];
            end
        end
      
 
        function output = getDrawableOutput(obj)
            % Rename default detection output
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Sub-resolution objects';
        end
        
        function hfigure = resultDisplay(obj,fig,procID)
            % Display the output of the process
              
            % Check for movie output before loading the GUI
            iChan = find(obj.checkChannelOutput,1);         
            if isempty(iChan)
                warndlg('The current step does not have any output yet.','No Output','modal');
                return
            end
            
            % Make sure detection output is valid
            movieInfo=obj.loadChannelOutput(iChan,'output','movieInfo');
            firstframe=find(arrayfun(@(x) ~isempty(x.amp),movieInfo),1);
            if isempty(firstframe)
                warndlg('The detection result is empty. There is nothing to visualize.','Empty Output','modal');
                return
            end
            
            hfigure = detectionVisualGUI('mainFig', fig, procID);
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'Gaussian Mixture-Model Fitting';
        end
        function h = GUI()
            h = @subResolutionProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % moviePara  
            funParams.ChannelIndex =1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'GaussianMixtureModels'];
            funParams.firstImageNum = 1;
            funParams.lastImageNum = owner.nFrames_;
            
            % detectionParam
            if ~isempty(owner.channels_(1).psfSigma_)
                funParams.detectionParam.psfSigma = owner.channels_.psfSigma_;
            else
                funParams.detectionParam.psfSigma=[];
            end
            if ~isempty(owner.camBitdepth_)
                funParams.detectionParam.bitDepth = owner.camBitdepth_;
            else
                funParams.detectionParam.bitDepth = [];
            end
            funParams.detectionParam.alphaLocMax = .05;
            funParams.detectionParam.integWindow = 0;
            funParams.detectionParam.doMMF = 0;
            funParams.detectionParam.testAlpha = struct('alphaR', .05,'alphaA', .05, 'alphaD', .05,'alphaF',0);
            funParams.detectionParam.numSigmaIter = 0;
            funParams.detectionParam.visual = 0;
            funParams.detectionParam.background = [];
            
        end

    end
    
end