classdef PointSourceDetectionProcess < DetectionProcess
    % A concrete class of a point source detection process
    
    methods (Access = public)
        function obj = PointSourceDetectionProcess(owner, outputDir, funParams )
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = PointSourceDetectionProcess.getName;
            super_args{3} = @detectMoviePointSources;
            
            if nargin < 3 || isempty(funParams)  % Default funParams
                if nargin <2, outputDir = owner.outputDirectory_; end
                funParams=PointSourceDetectionProcess.getDefaultParams(owner,outputDir);
            end
            
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
        end
        
        function movieInfo = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            ip =inputParser;
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output', 'movieInfo',@(x) strcmp(x, 'movieInfo'));
            ip.parse(iChan,varargin{:})
            
            % Data loading
            s = load(obj.outFilePaths_{1,iChan}, 'movieInfo');
            movieInfo = s.movieInfo(ip.Results.iFrame);
        end
        
        function output = getDrawableOutput(obj)
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Point sources';
            output(1).formatData=@PointSourceDetectionProcess.formatOutput;
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Point source detection';
        end
        function h = GUI()
            h = @pointSourceDetectionProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex=1;
            funParams.MaskChannelIndex = []; %1:numel(owner.channels_);
            funParams.MaskProcessIndex = [];            
            funParams.OutputDirectory = [outputDir  filesep 'point_sources'];
            funParams.alpha=.05;
            funParams.maskRadius=40;
            funParams.Mode = 'xyAc';
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RedundancyRadius = .25;
            funParams.UseIntersection = true;
            funParams.filterSigma = 1.2;
%             funParams.filterSigma = zeros(1, numel(owner.channels_));
%             hasPSFSigma = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
%             funParams.filterSigma(hasPSFSigma) = [owner.channels_(hasPSFSigma).psfSigma_];            
        end
        
        function positions = formatOutput(pstruct)
            positions = formatOutput@DetectionProcess(pstruct);
            %positions = positions(pstruct.isPSF, :);
        end
    end    
end