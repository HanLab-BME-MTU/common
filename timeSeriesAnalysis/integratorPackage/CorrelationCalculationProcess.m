classdef CorrelationCalculationProcess < CorrelationProcess
    % Process which
    %
    % Sebastien Besson, Oct 2011

    methods (Access = public)
        
        function obj = CorrelationCalculationProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                
                super_args{1} = owner;
                
                super_args{2} = CorrelationCalculationProcess.getName;
                super_args{3} = @calculateMovieCorrelation;
                
                if nargin < 3 || isempty(funParams)
                    funParams.OutputDirectory = [outputDir  filesep 'correlation'];
                    if isa(owner,'MovieList')
                        funParams.MovieIndex=1:numel(owner.movies_);
                    end
                    funParams.ProcessName=CorrelationProcess.getCorrelationProcesses;
                end
                
                super_args{4} = funParams;
                
            end
            
            obj = obj@CorrelationProcess(super_args{:});
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
        function procNames = getCorrelationProcesses()
            procNames = {'WindowSamplingProcess';
                'ProtrusionSamplingProcess'};
        end
    end
    
end


function data =formatCorrelationData(data)
data.X=data.lags;
data.Z=data.corrFun;
data=rmfield(data,{'lags','corrFun'});
end

