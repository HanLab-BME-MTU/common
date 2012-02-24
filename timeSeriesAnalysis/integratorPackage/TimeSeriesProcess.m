classdef TimeSeriesProcess < Process
    % Generic class to use for the time-series analysis
    %
    % Sebastien Besson, 7/2010 (last modified Mar, 2012)
    
    methods
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
        
        function input = getInput(obj,index)
            % Read process names from parameters
            procNames =obj.funParams_.ProcessName;
            nProc = numel(procNames);
            
            % Initialize process status
            procIndex = zeros(nProc,1);
            outputList = cell(nProc,1);
            isMovieProc = false(nProc,1);
            procOutput = cell(nProc,1);
            
            if isa(obj.owner_,'MovieList');
                movie=obj.owner_.movies_{1}; % Quick fix for movie lists
            else
                movie=obj.owner_;
            end
                
            % For each input process check the output validity
            for i=1:nProc
                procIndex(i) =movie.getProcessIndex(procNames{i},1);
                proc =movie.processes_{procIndex(i)};
                outputList{i} = proc.getDrawableOutput;
                isMovieProc(i) = strcmp('movieGraph',outputList{i}(1).type);
                procOutput{i} = proc.checkChannelOutput;
                assert(any(procOutput{i}),[proc.getName ' has no valid output !' ...
                    'Please apply ' proc.getName ' before running correlation!']);             
            end
            
            % Push all input into a structre
            nInput = sum(cellfun(@(x)sum(x(:)),procOutput));
            if nInput==0, input=[]; return; end
            input(nInput,1)=struct(); % Initialize time-series input structure
            iInput=0;
            for iProc=1:nProc
                for iOutput = 1:size(procOutput{iProc},1)
                    if isMovieProc(iProc)
                        % Add processIndex and output variable/name
                        iInput=iInput+1;
                        input(iInput).processIndex = procIndex(iProc);
                        input(iInput).var = outputList{iProc}(iOutput).var;
                        input(iInput).channelIndex = [];
                        input(iInput).name = regexprep(outputList{iProc}(iOutput).name,' map','');
                    else
                        % Loop over channels with valid output
                        for iChan=find(procOutput{iProc}(iOutput,:))
                            iInput=iInput+1;
                            input(iInput).processIndex = procIndex(iProc);
                            input(iInput).var = outputList{iProc}(iOutput).var;
                            input(iInput).outputIndex = iOutput;
                            input(iInput).channelIndex = iChan;
                            input(iInput).name = [regexprep(outputList{iProc}(iOutput).name,' map','') ' channel '...
                                num2str(iChan)];
                        end
                    end
                end
            end
            if nargin>1
                assert(all(ismember(index,1:numel(input))));
                input=input(index);
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
        
       
    end
    methods (Static)
        function procNames = getSamplingProcesses()
            procNames = {
                'ProtrusionSamplingProcess';
                'WindowSamplingProcess';};
        end
    end
end


