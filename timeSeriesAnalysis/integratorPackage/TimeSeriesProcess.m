classdef TimeSeriesProcess < Process
    % Process
    %
    % Sebastien Besson
    % 7/2010
    %
    methods (Access = public)
        
        function obj = TimeSeriesProcess(owner,name,funName,funParams)                         
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;                
            end
            
            obj = obj@Process(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;                              
            end
            if nargin > 3
               obj.funParams_ = funParams;              
            end
            
        end
        
        function input = getInput(obj,varargin)
            procNames =obj.funParams_.ProcessName;
            nProc = numel(procNames);
            
            procIndex = zeros(nProc,1);
            channelIndex=cell(nProc,1);
            outputList = cell(nProc,1);
            
            if isa(obj.owner_,'MovieList');
                movie=obj.owner_.movies_{1}; % Quick fix for movie lists
            else
                movie=obj.owner_;
            end
                
            % For each process check the channel/movie output
            for i=1:nProc
                procIndex(i) =movie.getProcessIndex(procNames{i},1);
                proc =movie.processes_{procIndex(i)};
                outputList{i} = proc.getDrawableOutput;
                isMovieProc = strcmp('movieGraph',outputList{i}(1).type);
                if isMovieProc
                    if ~proc.checkChannelOutput
                        error([proc.getName ' has no valid output !' ...
                            'Please apply ' proc.getName ' before running correlation!']);
                    end
                    channelIndex{i}=[];                    
                else
                    if ~any(proc.checkChannelOutput)
                        error([proc.getName ' has no valid output !' ...
                            'Please apply ' proc.getName ' before running correlation!']);
                    end
                    channelIndex{i}=find(proc.checkChannelOutput());  
                end
            end
            
            % Put all input in a structre
            procInNr = cellfun(@numel,channelIndex)+cellfun(@isempty,channelIndex);
            if isempty(procInNr), input=[]; return; end
            input(sum(procInNr))=struct();
            for i=1:nProc
                for j=1:procInNr(i)
                    inIdx=sum(procInNr(1:i-1))+j;
                    input(inIdx).processIndex = procIndex(i);
                    input(inIdx).var = outputList{i}.var;
                    if isempty(channelIndex{i})
                        input(inIdx).channelIndex = [];
                        input(inIdx).name = regexprep(outputList{i}.name,' map','');
                    else
                        input(inIdx).channelIndex = channelIndex{i}(j);
                        input(inIdx).name = [regexprep(outputList{i}.name,' map','') ' channel '...
                            num2str(channelIndex{i}(j))];
                    end
                end
            end
            if nargin>1
                assert(all(ismember(varargin{1},1:numel(input))));
                input=input(varargin{1});
            end
        end  
        
        function status = checkOutput(obj,varargin)
            % Input check
            input=obj.getInput;
            nInput=numel(input);
            ip =inputParser;
            ip.addOptional('iInput1',1:nInput,@(x) all(ismember(x,1:nInput)));
            ip.addOptional('iInput2',1:nInput,@(x) all(ismember(x,1:nInput)));
            ip.parse(varargin{:});
            iInput1=ip.Results.iInput1;
            iInput2=ip.Results.iInput2;
            
            %Makes sure there's at least one output file per channel
            status =  arrayfun(@(i,j) exist(obj.outFilePaths_{i,j},'file'),iInput1,iInput2);
    
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
            
            data=obj.loadOutput(i,j,'output',ip.Results.output);
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

    end
    methods (Static)
        function procNames = getSamplingProcesses()
            procNames = {
                'ProtrusionSamplingProcess';
                'WindowSamplingProcess';};
        end
    end
end


