classdef CorrelationCalculationProcess < TimeSeriesProcess
    % A concrete process for calculating correlation of sampled processes
    %
    % Sebastien Besson, Oct 2011 (last modified Dec 2011)
    
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
            else
                crossOutputList={'lags','ccf','ccfBounds',...
                    'bootstrapCcf','bootstrapCcfBounds'};
                partialOutputList = {'lags'};
            end
            outputList=unique(horzcat(crossOutputList,partialOutputList));
            
            % Parser optional output param/value
            additionalOutput=horzcat({''},{'cross'},{'partial'});
            allOutput = horzcat(additionalOutput,outputList);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,allOutput)));
            ip.parse(i,j,varargin{:});
            if strcmpi(ip.Results.output,'');
                output=outputList;
            elseif strcmpi(ip.Results.output,'cross');
                output=crossOutputList;
            elseif strcmpi(ip.Results.output,'partial');
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
            output(1).name='Cross correlation';
            output(1).var='cross';
            output(1).formatData=@formatCorrelation;
            output(1).type='correlationGraph';
            output(1).defaultDisplayMethod = @(i,j)ErrorBarGraphDisplay('XLabel','Lags (s)',...
                'YLabel',obj.getDrawableOutputName(i,j,'cross'),'Input1',obj.getInput(i).name);
            output(2).name='Partial correlation';
            output(2).var='partial';
            output(2).formatData=@formatCorrelation;
            output(2).type='correlationGraph';
            output(2).defaultDisplayMethod = @(i,j)ErrorBarGraphDisplay('XLabel','Lags (s)',...
                'YLabel',obj.getDrawableOutputName(i,j,'partial'),'Input1',obj.getInput(i).name);          
%             if isa(obj.owner_,'MovieData'),
%                 output(2).name='Correlation function';
%                 output(2).var='raw';
%                 output(2).formatData=@formatCorrelationData;
%                 output(2).type='correlationGraph';
%                 output(2).defaultDisplayMethod = @(i,j) MeshDisplay('XLabel','Lags (s)',...
%                     'YLabel','Window number','ZLabel',obj.getDrawableOutputName(i,j));
%             end
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
            funParams.nBoot=1e4;
            funParams.alpha=.01;
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
