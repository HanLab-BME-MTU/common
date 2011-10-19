classdef MaskIntersectionProcess < MaskProcessingProcess
    % A concrete process to create mask interesections
    
    
    methods(Access = public)        
        function obj = MaskIntersectionProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = MaskIntersectionProcess.getName;
                super_args{3} = @intersectMovieMasks;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%                                            
                    funParams.ChannelIndex = 1:nChan; %Default is to transform masks for all channels
                    funParams.SegProcessIndex = []; %No default...
                    funParams.OutputDirectory =  [outputDir filesep 'intersected_masks'];      
                    funParams.BatchMode = false;                                              
                    
                end
                
                super_args{4} = funParams;                    
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
            
        end  
    
        
    end
    methods(Static)
        function name =getName()
            name = 'Mask Intersection';
        end
    end
end
    