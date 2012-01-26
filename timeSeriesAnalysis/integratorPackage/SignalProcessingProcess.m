classdef SignalProcessingProcess < TimeSeriesProcess
    % A concrete process for calculating correlation of sampled processes
    %
    % Sebastien Besson, Oct 2011 (last modified Dec 2011)
    
    methods (Access = public)
        
        function obj = SignalProcessingProcess(owner,varargin)
            
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
                super_args{2} = SignalProcessingProcess.getName;
                super_args{3} = @processMovieSignal;
                if isempty(funParams)
                    funParams=SignalProcessingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@TimeSeriesProcess(super_args{:});
        end
        

        function varargout = loadOutput(obj,i,j,varargin)

            % Check input
            ip=inputParser;
            ip.addRequired('i',@isscalar);
            ip.addRequired('j',@isscalar);
            ip.parse(i,j);
            
            % Create specific output list for auto or cross correlation
            if i==j
                crossOutputList = {'lags','acf','acfBounds',...
                    'bootstrapAcf','bootstrapAcfBounds'};
                partialOutputList = {'lags','pacf','pacfBounds',...
                    'bootstrapPacf','bootstrapPacfBounds'};
                coherenceOutputList = {'avgCoh','cohCI','f'};
            else
                crossOutputList={'lags','ccf','ccfBounds',...
                    'bootstrapCcf','bootstrapCcfBounds'};
                partialOutputList = {'lags'};
            end
            outputList=unique(horzcat(crossOutputList,partialOutputList,coherenceOutputList));
            
            % Parser optional output param/value
            additionalOutput=horzcat({''},{'crossCorrelation'},{'partialCorrelation'},{'coherence'});
            allOutput = horzcat(additionalOutput,outputList);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,allOutput)));
            ip.parse(i,j,varargin{:});
            if strcmpi(ip.Results.output,'');
                output=outputList;
            elseif strcmpi(ip.Results.output,'crossCorrelation');
                output=crossOutputList;
            elseif strcmpi(ip.Results.output,'partial');
                output=partialOutputList;
            elseif strcmpi(ip.Results.output,'coherence');
                output=partialOutputList;  
            else
                output=ip.Results.output;
                if ischar(output), output={output}; end
            end

            % Load the data
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

            processingTools = obj.funParams_.processingTools;
            nOutput= numel(processingTools);
            output(nOutput,1)=struct();
            for i=1:numel(obj.funParams_.processingTools)
                output(i).name = tools(i).name;
                output(i).var = tools(i).var;
                output(i).formatData=tools(i).formatData;
                output(i).type = 'signalGraph';
                output(i).defaultDisplayMethod = tools(i).defaultDisplayMethod;
            end
        end
        
        function [label,title] = getDrawableOutputName(obj,i,j,output)
            switch output
                case 'cross'
                    if i==j,
                        label='Autocorrelation function';
                        title = [obj.getInput(i).name ' autocorrelation'];
                    else
                        label='Cross-correlation function';
                        title = [obj.getInput(i).name '/' obj.getInput(j).name ' cross-correlation'];
                    end
                case 'partial'
                    if i==j,
                        label='Partial autocorrelation function';
                        title = [obj.getInput(i).name ' partial autocorrelation'];
                    else
                        label='Partial cross-correlation function';
                        title = [obj.getInput(i).name '/' obj.getInput(j).name ' partial cross-correlation'];
                    end
                case 'coherence'
                    if i==j,
                        label='Power spectrum (db/Hz)';
                        title = [obj.getInput(i).name ' power spectrum'];
                    else
                        label='Coherence';
                        title = [obj.getInput(i).name '/' obj.getInput(j).name ' coherence'];
                    end
            end
        end
    end
    
    methods (Static)
        function name =getName()
            name = 'Correlation Calculation';
        end
        function h =GUI()
            h = @correlationCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'correlation'];
            funParams.ProcessName=TimeSeriesProcess.getTimeSeriesProcesses;
            if isa(owner,'MovieList'),
                funParams.MovieIndex=1:numel(owner.movies_);
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1,false)},...
                    owner.movies_,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) true(x.nSliceMax_,1),winProc,'UniformOutput',false);
            else
                winProc =owner.processes_{owner.getProcessIndex('WindowingProcess',1,false)};
                funParams.BandMin=1;
                funParams.BandMax=winProc.nBandMax_;
                funParams.SliceIndex=true(winProc.nSliceMax_,1);
            end
            
            tools = SignalProcessingProcess.getProcessingTools;
            funParams.processingTools = tools(1);
        end
        function tools = getProcessingTools()
            % Correlation function
            tools(1).name = 'Correlation';
            tools(1).GUI = @correlationSettingsGUI;
            tools(1).function = @measureDataCorrelation;
            tools(1).output = 'crossCorrelation';
            tools(1).outputList = 'crossCorrelation';
            tools(1).formatData = @formatCorrelation;
            tools(1).parameters.nBoot = 1e4;
            tools(1).parameters.alpha = .01;
            tools(1).defaultDisplayMethod = @(i,j)ErrorBarGraphDisplay('XLabel','Lags (s)',...
                    'YLabel',obj.getDrawableOutputName(i,j,'cross'),'Input1',obj.getInput(i).name);
            % Partial correlation function
            tools(2).name = 'Partial correlation';
            tools(2).GUI = @correlationSettingsGUI;
            tools(2).function = @measureDataPartialCorrelation;
            tools(2).formatData = @formatCorrelation;
            tools(2).parameters.nBoot = 1e4;
            tools(2).parameters.alpha = .01;
            % Coherence function
            tools(3).name = 'Coherence';
            tools(3).GUI = @coherenceSettingsGUI;
            tools(3).function = @measureDataCoherence;   
            tools(3).formatData = @formatCoherence;
            tools(3).parameters.nWin=8;
            tools(3).parameters.window='hamming';
            tools(3).parameters.noLap=.5;
            tools(3).parameters.nBoot=1e4;
            tools(3).parameters.alpha=.01;
            % Mutual information coefficient
            tools(4).name = 'Mutual information';
            tools(4).GUI = @micSettingsGUI;
            tools(4).function = @measureDataMutualInformation;   
            tools(4).parameters.nBoot=1e4;
        end
    end
    
end

function data =formatCorrelation(data)
data.X=data.lags;
if isfield(data,'bootstrapCcf')
    data.Y=data.bootstrapCcf;
    data.bounds=data.bootstrapCcfBounds;
elseif isfield(data,'bootstrapAcf')
    data.Y=data.bootstrapAcf;
    data.bounds=data.bootstrapAcfBounds;
else 
    data.Y=data.bootstrapPacf;
    data.bounds=data.bootstrapPacfBounds;
end
end

function data =formatCoherence(data)
data.X=data.f;
data.Y=data.avgCoh;
data.bounds=data.cohCI;
end
