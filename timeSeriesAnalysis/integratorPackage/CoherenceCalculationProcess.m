classdef CoherenceCalculationProcess < TimeSeriesProcess
    % A concrete process for calculating coherence of sampled processes
    %
    % Sebastien Besson, Dec 2011
    
    methods (Access = public)
        
        function obj = CoherenceCalculationProcess(owner,varargin)
            
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
                super_args{2} = CoherenceCalculationProcess.getName;
                super_args{3} = @calculateMovieCoherence;
                if isempty(funParams)
                    funParams=CoherenceCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@TimeSeriesProcess(super_args{:});
        end
        
        
        function varargout = loadChannelOutput(obj,i,j,varargin)
            
            % Check input
            ip=inputParser;
            ip.addRequired('i',@isscalar);
            ip.addRequired('j',@isscalar);
            ip.parse(i,j);
            
            
            % Check input
            outputList={'avgCoh','cohCI','f'};
            additionalOutput=horzcat({''});
            allOutput = horzcat(additionalOutput,outputList);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,allOutput)));
            ip.parse(i,j,varargin{:});
            if strcmpi(ip.Results.output,'');
                output=outputList;
            else
                output=ip.Results.output;
                if ischar(output), output={output}; end
            end
            
            s=load(obj.outFilePaths_{i,j},output{:});
            
            % Return structure if alll arguments are queried
            if any(ismember(ip.Results.output,additionalOutput)) || ...
                    isequal(ip.Results.output,outputList);
                varargout{1}=s;
            else
                for j=1:numel(output)
                    varargout{j} = s.(output{j});
                end
            end
            
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Coherence';
            output(1).var='';
            output(1).formatData=@formatCoherence;
            output(1).type='correlationGraph';
            output(1).defaultDisplayMethod = @(i,j)ErrorBarGraphDisplay('XLabel','Frequency (Hz)',...
                'YLabel',obj.getDrawableOutputName(i,j),'Input1',obj.getInput(i).name);
            
        end
        function [label,title] = getDrawableOutputName(obj,i,j,varargin)
            if i==j, 
                label='Power spectrum (db/Hz)'; 
                title = [obj.getInput(i).name ' power spectrum'];
            else
                label='Coherence'; 
                title = [obj.getInput(i).name '/' obj.getInput(j).name ' coherence'];
            end
        end
        
    end
    
    methods (Static)
        function name =getName()
            name = 'Coherence Calculation';
        end
        function h =GUI()
            h = @coherenceCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'coherence'];
            funParams.ProcessName=TimeSeriesProcess.getTimeSeriesProcesses;
            if isa(owner,'MovieList'),
                funParams.MovieIndex=1:numel(owner.movies_);
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1,false)},...
                    owner.movies_,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) ones(x.nSliceMax_,1),winProc,'UniformOutput',false);
            else
                winProc =owner.processes_{owner.getProcessIndex('WindowingProcess',1,false)};
                funParams.BandMin=1;
                funParams.BandMax=winProc.nBandMax_;
                funParams.SliceIndex=ones(winProc.nSliceMax_,1);
            end
            funParams.nWin=8;
            funParams.window='hamming';
            funParams.noLap=.5;
            funParams.nBoot=1e4;
            funParams.alpha=.01;
        end
    end
end

function data =formatCoherence(data)
data.X=data.f;
data.Y=data.avgCoh;
data.bounds=data.cohCI;
end
