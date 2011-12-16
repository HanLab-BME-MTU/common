classdef SignalPreprocessingProcess < TimeSeriesProcess
    % A concrete process for pre-processing time series
    %
    % Sebastien Besson, Oct 2011

    methods (Access = public)
        
        function obj = SignalPreprocessingProcess(owner,varargin)
            
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
                super_args{2} = SignalPreprocessingProcess.getName;
                super_args{3} = @preprocessMovieSignal;                
                if isempty(funParams)
                    funParams=SignalPreprocessingProcess.getDefaultParams(owner,outputDir);
                end                
                super_args{4} = funParams;                
            end
            
            obj = obj@TimeSeriesProcess(super_args{:});
        end
              
        
        function varargout = loadChannelOutput(obj,i,varargin)
            % Check input
            outputList={'data','range'};
            ip=inputParser;
            ip.addRequired('obj');
            ip.addRequired('i',@isscalar);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(obj,i,varargin{:});
            output=ip.Results.output;
            if ischar(output), output={output}; end
            
            s=load(obj.outFilePaths_{1,i},output{:});
            for j=1:numel(output)
                if strcmp(output{j},'')
                    varargout{j}=s;
                else
                    varargout{j} = s.(output{j});
                end
            end
        end
        
        function status = checkChannelOutput(obj,i)
            status = cellfun(@(x)exist(x,'file'),obj.outFilePaths_(1,i));
        end  
    end
    
    methods (Static)
        function name =getName()
            name = 'Signal Preprocessing';
        end
        function h =GUI()
            h = @signalPreprocessingProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            if isa(owner,'MovieList'), funParams.MovieIndex=1:numel(owner.movies_); end
            funParams.OutputDirectory = [outputDir  filesep 'preprocessedSignal'];
            funParams.ProcessName=SignalPreprocessingProcess.getTimeSeriesProcesses;
            funParams.kSigma=5;
        end
    end
end
