classdef DummySegmentationProcess < SegmentationProcess
    % A concrete class for importing external mask output
    
    methods(Access = public)
        
        function obj = DummySegmentationProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            
            % Constructor of the DummyDetectionProcess
            super_args{1} = owner;
            super_args{2} = DummySegmentationProcess.getName();
            super_args{3} = @runDummySegmentation;
            if isempty(ip.Results.funParams)
                super_args{4} = DummySegmentationProcess.getDefaultParams(...
                    owner, ip.Results.outputDir);
            else
                super_args{4} = ip.Results.funParams;
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'Dummy segmentation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'dummySegmentation'];
            funParams.InputData = cell(numel(owner.channels_), 1);
        end
    end
end
