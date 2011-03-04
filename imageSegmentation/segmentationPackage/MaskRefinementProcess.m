classdef MaskRefinementProcess < MaskProcessingProcess
    %Class definition for post processing using refineMovieMasks.m
    
    
    methods(Access = public)        
        function obj = MaskRefinementProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = 'Mask Refinement';
                super_args{3} = @refineMovieMasks;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%                                            
                    funParams.ChannelIndex = 1:nChan; %Default is to attempt to refine masks for all channels
                    funParams.SegProcessIndex = []; %No default.
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'refined_masks'];      
                    funParams.MaskCleanUp = true;
                    funParams.MinimumSize = 10;
                    funParams.ClosureRadius = 3;
                    funParams.ObjectNumber = 1; %Default is to keep only 1 object per mask
                    funParams.FillHoles = true;
                    funParams.EdgeRefinement = false; %This off by default because it sort of sucks, and is slow.
                    funParams.MaxEdgeAdjust = []; %Use refineMaskEdges.m function defaults for these settings
                    funParams.MaxEdgeGap = [];    
                    funParams.PreEdgeGrow = [];
                    funParams.BatchMode = false;                                              
                end
                %Make sure the input parameters are legit??
                super_args{4} = funParams;                    
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
            
        end                  
        
        
        
    end
end
    