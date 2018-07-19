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
              
        
        function varargout = loadOutput(obj,i,varargin)
            % Check input
            outputList={'data','range','energyData','energy'};
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
        
        function h=draw(obj,iInput,varargin)
            
            % Check input
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('iInput',@isscalar);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(iInput,varargin{:})
			
            % Load data
            data=obj.loadOutput(iInput,'output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
           
            try
                assert(~isempty(obj.displayMethod_{iInput}));
            catch ME %#ok<NASGU>
                obj.displayMethod_{iInput}=outputList(iOutput).defaultDisplayMethod(iInput);
            end
            
            % Delegate to the corresponding method
            tag = ['process' num2str(obj.getIndex) '_input' num2str(iInput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            input=obj.getInput;
            procArgs={'Input1',input(1).name};
            h=obj.displayMethod_{iInput}.draw(data,tag,drawArgs{:},procArgs{:});
        end
        
        function status = checkOutput(obj,varargin)
            % Input check
            input=obj.getInput;
            nInput=numel(input);
            ip =inputParser;
            ip.addOptional('iInput',1:nInput,@(x) all(ismember(x,1:nInput)));
            ip.parse(varargin{:});
            iInput=ip.Results.iInput;
            
            if isempty(obj.outFilePaths_), 
                status =false(size(iInput));
            else
                status = cellfun(@(x) logical(exist(x,'file')),obj.outFilePaths_(iInput));
            end
        end  
        
        function output = getDrawableOutput(obj)
            output.name = 'Energy';
            output.var = 'energy';
            output.formatData=@formatEnergy;
            output.type = 'signalGraph';
            output.defaultDisplayMethod = @(i)ErrorBarGraphDisplay('XLabel','Band',...
                'YLabel','Mean energy','Input1',obj.getInput(i).name); 
        end
        
    end
    
    methods (Static)
        function name =getName()
            name = 'Signal Preprocessing';
        end
        function h =GUI()
            h = @signalPreprocessingProcessGUI;
        end
        function trends= getTrends()
            names = {'None','Mean','Linear','Deterministic'};
            type = {-1,0,1,2};
            trends=cell2struct(vertcat(names,type),{'name','type'});
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            if isa(owner,'MovieList'), 
                signal=owner.getMovies{1}.getSampledOutput;
                funParams.MovieIndex=1:numel(owner.getMovies); 
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1,false)},...
                    owner.getMovies,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) true(x.nSliceMax_,1),winProc,'UniformOutput',false);
   
            else
                signal=owner.getSampledOutput;
                winProcIndex = owner.getProcessIndex('WindowingProcess',1,false);
                assert(~isempty(winProcIndex),'lccb:getDefaultParams:noWindowing','Movie is not windowed yet');
                winProc =owner.processes_{winProcIndex};
                funParams.BandMin=1;
                funParams.BandMax=winProc.nBandMax_;
                funParams.SliceIndex=true(winProc.nSliceMax_,1);
            end
            nSignal=numel(signal);
            
            funParams.OutputDirectory = [outputDir filesep 'preprocessedSignal'];       
            funParams.kSigma=5*ones(nSignal,1);
            funParams.trendType=1*ones(nSignal,1);
            funParams.nBoot=1e3;
            funParams.alpha=.05;
        end
        
    end
end

function data =formatEnergy(energy)
data.X=1:size(energy,1);
data.Y=energy(:,1);
data.bounds=energy(:,2:3)';
end

