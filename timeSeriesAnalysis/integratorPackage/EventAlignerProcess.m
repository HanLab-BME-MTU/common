classdef EventAlignerProcess < Process
    % Process
    %
    % Sebastien Besson Oct 2011
    methods (Access = public)
        
        function obj = EventAlignerProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                
                super_args{1} = owner;
                
                super_args{2} = CorrelationBootstrappingProcess.getName;
                super_args{3} = @bootstrapMoviesCorrelation;
                
                if nargin < 3 || isempty(funParams)    
                    funParams.OutputDirectory = [outputDir  filesep 'eventAlignmnet'];
 
                end
                
                super_args{4} = funParams;
                
            end
            
            obj = obj@Process(super_args{:});
        end

    end
    
    methods (Static)
        function name =getName()
            name = 'Event Aligner';
        end


    end
    
end

