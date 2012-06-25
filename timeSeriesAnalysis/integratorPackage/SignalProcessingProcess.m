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
        
        function status = checkOutput(obj,varargin)
            % Input check
            input=obj.getInput;
            nInput=numel(input);
            nTools = numel(obj.funParams_.tools);
            ip =inputParser;
            ip.addOptional('iInput1',1:nInput,@(x) all(ismember(x,1:nInput)));
            ip.addOptional('iInput2',1:nInput,@(x) all(ismember(x,1:nInput)));
            ip.addOptional('iOutput',1:nTools,@(x) all(ismember(x,1:nTools)));
            ip.parse(varargin{:});
            iInput1=ip.Results.iInput1;
            iInput2=ip.Results.iInput2;
            iOutput = ip.Results.iOutput;
            
            %Makes sure there's at least one output file per channel
            status = cellfun(@(x) logical(exist(x,'file')),obj.outFilePaths_(iInput1,iInput2,iOutput));
        end

        function varargout = loadOutput(obj,iInput1,iInput2,iOutput)

            % Check input
            input=obj.getInput;
            nInput=numel(input);
            tools=obj.funParams_.tools;
            nTools = numel(tools);
            ip=inputParser;
            ip.addRequired('iInput1',@(x) isscalar(x) && ismember(x,1:nInput));
            ip.addRequired('iInput2',@(x) isscalar(x) && ismember(x,1:nInput));
            ip.addRequired('iOutput',@(x) isscalar(x) && ismember(x,1:nTools));
            ip.parse(iInput1,iInput2,iOutput);
            
            % Create specific output list for auto or cross correlation
            allTools = obj.getTools;
            selectedTool = allTools(tools(iOutput).type);            
            
            % Load the data
            vars = whos('-file',obj.outFilePaths_{iInput1,iInput2,iOutput});
            outputList = intersect({vars.name},selectedTool.outputList);
            varargout{1}=load(obj.outFilePaths_{iInput1,iInput2,iOutput},outputList{:});
        end
        
        function h=draw(obj,iInput1,iInput2,iOutput,varargin)
            
            % Check input
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('iInput1',@isscalar);
            ip.addRequired('iInput2',@isscalar);
            ip.addRequired('iOutput',@isscalar);
            ip.KeepUnmatched = true;
            ip.parse(iInput1,iInput2,iOutput,varargin{:})

            data=obj.loadOutput(iInput1,iInput2,iOutput);
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            
            try
                assert(~isempty(obj.displayMethod_{iInput1,iInput2,iOutput}));
            catch ME %#ok<NASGU>
                obj.displayMethod_{iInput1,iInput2,iOutput}=outputList(iOutput).defaultDisplayMethod(iInput1,iInput2);
            end
            
            % Delegate to the corresponding method
            tag = ['process' num2str(obj.getIndex) '_input' num2str(iInput1) '_input' num2str(iInput2)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            input=obj.getInput;
            procArgs={'Input1',input(1).name,'Input2',input(2).name};
            h=obj.displayMethod_{iInput1,iInput2,iOutput}.draw(data,tag,drawArgs{:},procArgs{:});
        end

        
        function output = getDrawableOutput(obj)

            allTools = obj.getTools;
            tools = obj.funParams_.tools;
            nOutput= numel(tools);
            output(nOutput,1)=struct();
            
            for i=1:numel(tools)
                selectedTool = allTools(tools(i).type);
                output(i).name = selectedTool.name;
                output(i).var = selectedTool.name;
                output(i).formatData=selectedTool.formatData;
                output(i).type = 'signalGraph';
                output(i).defaultDisplayMethod = @(i1,i2)ErrorBarGraphDisplay('XLabel',selectedTool.xlabel,...
                'YLabel',selectedTool.ylabel(i1,i2),'Input1',obj.getInput(i1).name);
            end
        end
        
        function title=getOutputTitle(obj,iInput1,iInput2,iOutput)
            toolType =obj.funParams_.tools(iOutput).type;
            if iInput1==iInput2,
                switch toolType
                    case 1, name = 'autocorrelation';
                    case 2, name = 'partial autocorrelation';
                    case 3, name = 'power spectrum';
                    otherwise , error('Undefined output type');
                end
                title = [obj.getInput(iInput1).name ' ' name];
            else
                switch toolType
                    case 1, name = 'cross-correlation';
                    case 2, name = 'partial cross-correlation';
                    case 3, name = 'coherence';
                    otherwise , error('Undefined output type');
                end
                title = [obj.getInput(iInput1).name '/' obj.getInput(iInput2).name ' ' name];
            end
        end
        
        
        function addTool(obj,tools)
            ip=inputParser;
            ip.addRequired('tools',@isstruct)
            ip.parse(tools);
            
            fields=fieldnames(tools);
            obj.funParams_.tools(end+1)=rmfield(tools,fields(~ismember(fields,{'type','parameters'})));
        end
        
        function removeTool(obj,i)
            obj.funParams_.tools(i)=[];
        end
    end
    
    methods (Static)
        function name =getName()
            name = 'Signal Processing';
        end
        function h =GUI()
            h = @signalProcessingProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'signalProcessing'];
            if isa(owner,'MovieList'),
                funParams.MovieIndex=1:numel(owner.getMovies);
                signal=owner.getMovies{1}.getSampledOutput;  
            else
                signal=owner.getSampledOutput;
            end            
            funParams.SignalIndex=1:numel(signal);
            tools = SignalProcessingProcess.getTools(1);
            fields=fieldnames(tools);
            funParams.tools=rmfield(tools,fields(~ismember(fields,{'type','parameters'})));
        end
        function tools = getTools(index)
            % Correlation function
            tools(1).type = 1;
            tools(1).name = 'Correlation';            
            tools(1).GUI = @correlationSettingsGUI;
            tools(1).function = @computeSignalCorrelation;
            tools(1).outputList = {'lags','acf','acfBounds','bootstrapAcf',...
                'bootstrapAcfBounds','ccf','ccfBounds','bootstrapCcf',...
                'bootstrapCcfBounds'};
            tools(1).xlabel = 'Lags (s)';
            tools(1).ylabel = @(i,j) getCorrelationLabel(i,j);
            tools(1).title = @(i,j) getCorrelationTitle(i,j);
            tools(1).formatData = @formatCorrelation;
            tools(1).parameters = struct('nBoot',1e4,'alpha',.01);
            % Partial correlation function
            tools(2).type = 2;
            tools(2).name = 'Partial correlation';
            tools(2).outputList = {'lags','pacf','pacfBounds','bootstrapPacf',...
                'bootstrapPacfBounds'};
            tools(2).GUI = @correlationSettingsGUI;
            tools(2).function = @computeSignalPartialCorrelation;
            tools(2).formatData = @formatCorrelation;
            tools(2).xlabel = 'Lags (s)';
            tools(2).ylabel = @(i,j) getPartialCorrelationLabel(i,j);
            tools(2).title = @(i,j) getPartialCorrelationTitle(i,j);
            tools(2).parameters = struct('nBoot',1e4,'alpha',.01);
            % Coherence function
            tools(3).type = 3;
            tools(3).name = 'Coherence';
            tools(3).output = 'coherence';
            tools(3).outputList = {'f','avgCoh','cohCI'};
            tools(3).GUI = @coherenceSettingsGUI;
            tools(3).function = @computeSignalCoherence;   
            tools(3).formatData = @formatCoherence;
            tools(3).parameters = struct('nBoot',1e4,'alpha',.01,...
                'window','hamming','noLap',.5,'nWin',8);
            tools(3).xlabel = 'Frequency (Hz)';
            tools(3).ylabel = @(i,j) getCoherenceLabel(i,j);
            tools(3).title = @(i,j) getCoherenceTitle(i,j);
            
            if nargin>0, tools=tools(index); end
        end
        
        function windows = getCoherenceWindows()
            windows={'bartlett';
                'blackman';
                'blackmanharris';
                'bohman';
                'bohmanwin';
                'chebwin';
                'flattopwin';
                'gausswin';
                'hamming';
                'hann';
                'kaiser';
                'nuttallwin';
                'parzenwin';
                'rectwin';
                'taylorwin';
                'triang';
                'tukeywin'};
        end
    end
end

function label =getCorrelationLabel(i,j)
if i==j,
    label='Autocorrelation function';
else
    label='Cross-correlation function';
end
end

function label =getPartialCorrelationLabel(i,j)
if i==j, label='Partial autocorrelation function'; else label='Partial cross-correlation function';end
end

function label =getCoherenceLabel(i,j)
if i==j, label='Power spectrum (dB/Hz)'; else label='Coherence'; end
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
