classdef EventAlignmentProcess < TimeSeriesProcess
    % A concrete process for aligning events of sampled processes
    %
    % Sebastien Besson, Dec 2011
    
    methods (Access = public)
        
        function obj = EventAlignmentProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = EventAlignmentProcess.getName;
                super_args{3} = @alignMovieEvents;
                if isempty(funParams)
                    funParams=EventAlignmentProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@TimeSeriesProcess(super_args{:});
        end
        
        
        function varargout = loadChannelOutput(obj,i,j,varargin)
            % Check input
            outputList={'','corrFun','bounds','lags'};
            ip=inputParser;
            ip.addRequired('obj');
            ip.addRequired('i',@isscalar);
            ip.addRequired('j',@isscalar);
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,i,j,varargin{:});
            output=ip.Results.output;
            if ischar(output), output={output}; end
            
            s=load(obj.outFilePaths_{i,j},output{:});
            for j=1:numel(output)
                if strcmp(output{j},'')
                    varargout{j}=s;
                else
                    varargout{j} = s.(output{j});
                end
            end
        end
        
        
        function output = getDrawableOutput(obj)
            output(1).name='Correlation';
            output(1).var={''};
            output(1).formatData=[];
            output(1).type='correlationGraph';
            output(1).defaultDisplayMethod = @ScalarMapDisplay;
        end
    end
    
    methods (Static)
        function name =getName()
            name = 'Event Alignment';
        end
        function h =GUI()
            h = @eventAlignmentProcessGUI;
        end
        
        function events = getEvents(processname)
            switch(processname)
                case('ProtrusionSamplingProcess')
                    events(1).name='Maximum protrusion velocity';
                    events(2).name='Maximum retraction velocity';
            end
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'events'];
            funParams.ProcessName=TimeSeriesProcess.getTimeSeriesProcesses;
            if isa(owner,'MovieList'),
                funParams.MovieIndex=1:numel(owner.movies_);
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1)},...
                    owner.movies_,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) ones(x.nSliceMax_,1),winProc,'UniformOutput',false);
            end
            funParams.AlignmentProcess = [];
            funParams.Event = [];
        end
    end
end