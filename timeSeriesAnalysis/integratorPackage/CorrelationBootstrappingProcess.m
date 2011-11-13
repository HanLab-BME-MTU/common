classdef CorrelationBootstrappingProcess < CorrelationProcess
    % A concrete process for calculating confidence interval of correlation
    
    % Sebastien Besson Oct 2011
    methods (Access = public)
        
        function obj = CorrelationBootstrappingProcess(owner,varargin)
            
            if nargin>0
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = CorrelationBootstrappingProcess.getName;
                super_args{3} = @bootstrapMoviesCorrelation;
                if isempty(funParams)
                    funParams=CorrelationBootstrappingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@CorrelationProcess(super_args{:});
        end
        
        
        
        function varargout = loadChannelOutput(obj,i,j,varargin)
            % Check input
            outputList={'','avgCorrFun','avgBounds','lags'};
            ip=inputParser;
            ip.addRequired('obj');
            ip.addRequired('i',@isscalar);
            ip.addRequired('j',@isscalar);
            ip.addParamValue('output',outputList{1},@(x) any(strcmp(x,outputList)));
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
        
        
        function status = checkChannelOutput(obj,i,j)
            
        end
        
    end
    
    methods (Static)
        function name =getName()
            name = 'Correlation Bootstrapping';
        end
        function h =GUI()
            h = @correlationBootstrappingProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieList'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'bootstrappedCorrelation'];
            funParams.MovieIndex=1:numel(owner.movies_);
            funParams.ProcessName=CorrelationProcess.getCorrelationProcesses;
            winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1)},...
                owner.movies_,'UniformOutput',false);
            funParams.BandMin=1;
            funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
            funParams.SliceIndex=cellfun(@(x) 1:x.nSliceMax_,winProc,'UniformOutput',false);
        end
    end
end