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
            if isa(obj.owner_,'MovieList');
                input=obj.owner_.getMovies{1}.getSampledOutput; % Quick fix for movie lists
            else
                input=obj.owner_.getSampledOutput;
            end
            if nargin>1, input=input(index); end
                
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


