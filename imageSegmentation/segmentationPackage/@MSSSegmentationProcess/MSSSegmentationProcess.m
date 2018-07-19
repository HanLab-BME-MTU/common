classdef MSSSegmentationProcess < SegmentationProcess
    %A concrete process for segmenting using multi-scale steerable filters
    
    % Sebastien Besson, Sep 2011 (last modified Nov 2011)
    
    methods
        function obj = MSSSegmentationProcess(owner,varargin)
            
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
                super_args{2} = MSSSegmentationProcess.getName;
                super_args{3} = @getMovieMasksMSS;
                if isempty(funParams)
                    funParams=MSSSegmentationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Multi-scale segmentation';
        end
        function h = GUI()
            h= @mssSegmentationProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'MSSMasks'];
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.Scales = [1 2 4]; %Default is no jump suppression
            funParams.FilterOrder = 3;
        end
    end
end