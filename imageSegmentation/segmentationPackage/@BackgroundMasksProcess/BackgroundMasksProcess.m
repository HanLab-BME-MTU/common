classdef BackgroundMasksProcess < MaskProcessingProcess
    % A concrete process for creating background masks
    
    methods(Access = public)
        
        function obj = BackgroundMasksProcess(owner,varargin)
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
                super_args{2} = BackgroundMasksProcess.getName;
                super_args{3} = @createMovieBackgroundMasks;                
                if isempty(funParams)
                  funParams=BackgroundMasksProcess.getDefaultParams(owner,outputDir); 
                end
                super_args{4} = funParams;
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
        end
        
    end
    methods(Static)
        function name =getName()
            name = 'Background Mask';
        end
        function h = GUI()
            h= @backgroundMasksProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_); %Default is to attempt to creat background masks for all channels
            funParams.SegProcessIndex = []; %No default...
            funParams.OutputDirectory = [outputDir  filesep 'background_masks'];
            funParams.GrowthRadius = 20;
            funParams.BatchMode = false;
        end
    end
end