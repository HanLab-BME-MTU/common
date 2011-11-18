classdef CorrelationCalculationProcess < TimeSeriesProcess
    % A concrete process for calculating correlation of sampled processes
    %
    % Sebastien Besson, Oct 2011

    methods (Access = public)
        
        function obj = CorrelationCalculationProcess(owner,varargin)
            
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
                super_args{2} = CorrelationCalculationProcess.getName;
                super_args{3} = @calculateMovieCorrelation;                
                if isempty(funParams)
                    funParams=CorrelationCalculationProcess.getDefaultParams(owner,outputDir);
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
            output(1).formatData=@formatCorrelationData;
            output(1).type='correlationGraph';
            output(1).defaultDisplayMethod = @CorrelationMeshDisplay;
        end
    end
    
    methods (Static)
        function name =getName()
            name = 'Correlation Calculation';
        end
        function h =GUI()
            h = @correlationCalculationProcessGUI;
        end
        function procNames = getTimeSeriesProcesses()
            procNames = {'WindowSamplingProcess';
                'ProtrusionSamplingProcess'};
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieList'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'correlation'];
            funParams.ProcessName=TimeSeriesProcess.getTimeSeriesProcesses;
            if isa(owner,'MovieList'), 
                funParams.MovieIndex=1:numel(owner.movies_);
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1)},...
                    owner.movies_,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) 1:x.nSliceMax_,winProc,'UniformOutput',false);
            end
        end
    end
end

function data =formatCorrelationData(data)
data.X=data.lags;
data.Z=data.corrFun;
data=rmfield(data,{'lags','corrFun'});
end

