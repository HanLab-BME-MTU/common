classdef CorrelationProcess < Process
    % Process
    %
    % Sebastien Besson
    % 7/2010
    %
    methods (Access = public)
        
        function obj = CorrelationProcess(owner,name,funName,funParams)                         
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
        
        function input = getInput(obj)
            procNames =obj.funParams_.ProcessName;
            nProc = numel(procNames);
            
            procIndex = zeros(nProc,1);
            channelIndex=cell(nProc,1);
            outputList = cell(nProc,1);

            for i=1:nProc
                procIndex(i) = obj.owner_.getProcessIndex(procNames{i},1);
                proc =obj.owner_.processes_{procIndex(i)};
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
            
            % List input
            procInNr = cellfun(@numel,channelIndex)+cellfun(@isempty,channelIndex);
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
        end
        
        function status = checkChannelOutput(obj,i,j)
            status = cellfun(@(x)exist(x,'file'),obj.outFilePaths_(i,j));
        end  
        
        function h=draw(obj,i,varargin)
            % Function to draw process output (template method)
            
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
            input=obj.getInput;
            if ~isempty(outputList(1).formatData),
                data=outputList(1).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{i,j}));
            catch ME
                obj.displayMethod_{i,j}=outputList(1).defaultDisplayMethod();
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_process' num2str(i) '_process' num2str(j)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            procArgs={'Input1',input(1).name,'Input2',input(2).name};
            h=obj.displayMethod_{i,j}.draw(data,tag,drawArgs{:},procArgs{:});
        end
        function status = sanityCheck(obj)
                    
        end
        
    end
    methods (Static)
        function procNames = getCorrelationProcesses()
            procNames = {'ProtrusionSamplingProcess';
                'WindowSamplingProcess';};
        end
    end
end


