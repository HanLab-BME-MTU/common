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
        
        function input = getInput(obj)
            procNames =obj.funParams_.ProcessName;
            nProc = numel(procNames);
            
            procIndex = zeros(nProc,1);
            channelIndex=cell(nProc,1);
            outputList = cell(nProc,1);
            
            if isa(obj.owner_,'MovieList');
                movie=obj.owner_.movies_{1};
            else
                movie=obj.owner_;
            end
                
            
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
            
            % List input
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
        end  
        
        function [procIndex,channelIndex] =input2procchannel(obj,iInput)
            input=obj.getInput;
            procIndex = input(iInput).processIndex;
            channelIndex = input(iInput).channelIndex;
        end
        

    end
    methods (Static)
        function procNames = getTimeSeriesProcesses()
            procNames = {'ProtrusionSamplingProcess';
                'WindowSamplingProcess';};
        end
    end
end


