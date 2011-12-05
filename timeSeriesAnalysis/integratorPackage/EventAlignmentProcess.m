classdef EventAlignmentProcess < TimeSeriesProcess
    % A concrete process for aligning events of sampled processes
    %
    % Sebastien Besson, Dec 2011
    
    methods (Access = public)
        
        function obj = EventAlignmentProcess(owner,varargin)
            
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
                super_args{2} = EventAlignmentProcess.getName;
                super_args{3} = @alignMovieEvents;
                if isempty(funParams)
                    funParams=EventAlignmentProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@TimeSeriesProcess(super_args{:});
        end
        
        function varargout = loadChannelOutput(obj,iInput,varargin)
                                   
            %Input check
            input=obj.getInput;
            ip=inputParser;
            ip.addRequired('iInput',@isscalar);
            ip.addParamValue('output',input(iInput).var,@ischar);
            ip.parse(iInput,varargin{:});
            output = ip.Results.output;
            
            % Load
            s=load(obj.outFilePaths_{1,iInput},output);
            varargout{1}=s.(output);
        end
        
        function h=draw(obj,iInput,varargin)
            
            % Check input
            ip=inputParser;
            ip.addRequired('iInput',@isscalar);
            outputList = obj.getDrawableOutput(iInput);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(iInput,varargin{:})

            data=obj.loadChannelOutput(iInput,'output',ip.Results.output);     
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            
            try
                assert(~isempty(obj.displayMethod_{iOutput,iInput}));
            catch ME %#ok<NASGU>
                obj.displayMethod_{iOutput,iInput}=outputList(iOutput).defaultDisplayMethod(iInput);
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_output' num2str(iOutput) '_input' num2str(iInput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput,iInput}.draw(data,tag,drawArgs{:});
         end
        
        
                
        function status = checkChannelOutput(obj,varargin)
            
            %Input check
            input=obj.getInput;
            nInput=numel(input);
            ip=inputParser;
            ip.addOptional('iInput',1:nInput,@(x) all(ismember(x,1:nInput)));
            ip.parse(varargin{:});
            iInput=ip.Results.iInput;         
             
            status =  logical(arrayfun(@(x) exist(obj.outFilePaths_{1,x},...
                'file'),iInput)); 
        end  

        
        function output = getDrawableOutput(obj,varargin)
           
            output(1).name='Aligned maps';
            if nargin>1
                iInput=varargin{1};
                input=obj.getInput;
                original_output=obj.owner_.processes_{input(iInput).processIndex}.getDrawableOutput;
                output(1).var=input(iInput).var;
                output(1).formatData=original_output.formatData;
                output(1).type=original_output.type;
                output(1).defaultDisplayMethod = @(x) ScalarMapDisplay('Colormap','jet',...
                'CLim',obj.getIntensityLimits(x),'Labels',{'Frame number','Window depth','Window number'});
            else
                output(1).var='';
                output(1).formatData=[];
                output(1).type='graph';
                output(1).defaultDisplayMethod = []; 
            end
         

        end
    end
    
    methods (Access=protected)
        function limits = getIntensityLimits(obj,iInput)
            data=obj.loadChannelOutput(iInput);
            limits=[min(data(:)) max(data(:))];
        end
    end
    methods (Static)
        function name =getName()
            name = 'Event Alignment';
        end
        function h =GUI()
            h = @eventAlignmentProcessGUI;
        end
        
        function events = getEvents(processname)
            switch(processname)
                case('ProtrusionSamplingProcess')
                    events(1).name='Maximum protrusion velocity';
                    events(1).func=@max;
                    events(2).name='Minimum retraction velocity';
                    events(2).func=@min;
            end
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'events'];
            funParams.ProcessName=TimeSeriesProcess.getTimeSeriesProcesses;
            if isa(owner,'MovieList'),
                funParams.MovieIndex=1:numel(owner.movies_);
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1)},...
                    owner.movies_,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) ones(x.nSliceMax_,1),winProc,'UniformOutput',false);
            end
        end
    end
end