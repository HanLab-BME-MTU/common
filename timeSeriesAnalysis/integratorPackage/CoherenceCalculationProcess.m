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
            outputList={'raw','bootstrap','P','f'};
            ip=inputParser;
            ip.addRequired('obj');
            ip.addRequired('i',@isscalar);
            ip.addRequired('j',@isscalar);
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,i,j,varargin{:});
            output=ip.Results.output;
            if ischar(output), output={output}; end
            
            if strcmp(output{:},'raw')
                s=load(obj.outFilePaths_{i,j},'P','f');
            elseif strcmp(output{:},'bootstrap')
                s=load(obj.outFilePaths_{i,j},'bootstrapP','f');
            else
                s=load(obj.outFilePaths_{i,j},output{:});
            end
            
            for j=1:numel(output)
                if ismember(output{j},{'raw','bootstrap'})
                    varargout{j}=s;
                else
                    varargout{j} = s.(output{j});
                end
            end
        end
        
        function h=draw(obj,i,varargin)
            
            % Check input
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('i',@isscalar);
            ip.addOptional('j',i,@isscalar);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,i,varargin{:})
            j=ip.Results.j;
            
            data=obj.loadChannelOutput(i,j,'output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            
            try
                assert(~isempty(obj.displayMethod_{iOutput,i,j}));
            catch ME %#ok<NASGU>
                obj.displayMethod_{iOutput,i,j}=outputList(iOutput).defaultDisplayMethod(i,j);
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_input' num2str(i) '_input' num2str(j)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            input=obj.getInput;
            procArgs={'Input1',input(1).name,'Input2',input(2).name};
            h=obj.displayMethod_{iOutput,i,j}.draw(data,tag,drawArgs{:},procArgs{:});
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Bootsrapped spectral density';
            output(1).var='bootstrap';
            output(1).formatData=@formatBootstrappedSpectralDensity;
            output(1).type='correlationGraph';
            output(1).defaultDisplayMethod = @CorrelationGraphDisplay;
            if isa(obj.owner_,'MovieData'),
                output(2).name='Spectral density';
                output(2).var='raw';
                output(2).formatData=@formatSpectralDensity;
                output(2).type='correlationGraph';
                output(2).defaultDisplayMethod = @(i,j) MeshDisplay('XLabel','Frequency (Hz)',...
                    'YLabel','Window number','ZLabel',obj.getDrawableOutputName(i,j));
            end
            
        end
        function [label,title] = getDrawableOutputName(obj,i,j)
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
            funParams.window=[];
        end
    end
end

function data =formatSpectralDensity(data)
data.X=data.f;
data.Z=10*log10(data.P);
data=rmfield(data,{'P','f'});
end

function data =formatBootstrappedSpectralDensity(data)
data.lags=squeeze(nanmean(data.lags,2));
data.avgCorrFun=data.bootstrapCorrFun;
data.steCorrFun=data.bootstrapSteCorrFun;
end
