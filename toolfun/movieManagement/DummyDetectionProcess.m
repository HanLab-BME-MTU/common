classdef DummyDetectionProcess < DetectionProcess
    % A concrete class for importing external detection output
    
    methods(Access = public)
        
        function obj = DummyDetectionProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            
            % Constructor of the DummyDetectionProcess
            super_args{1} = owner;
            super_args{2} = DummyDetectionProcess.getName();
            super_args{3} = @runDummyDetection;
            if isempty(ip.Results.funParams)
                super_args{4} = DummyDetectionProcess.getDefaultParams(...
                    owner, ip.Results.outputDir);
            else
                super_args{4} = ip.Results.funParams;
            end
            
            obj = obj@DetectionProcess(super_args{:});
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'Dummy detection';
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
            funParams.OutputDirectory = [outputDir  filesep 'dummyDetection'];
            funParams.InputData = cell(numel(owner.channels_), 1);
        end
    end
end
