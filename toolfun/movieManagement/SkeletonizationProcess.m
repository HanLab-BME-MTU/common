classdef SkeletonizationProcess < ImageProcessingProcess
% A process class for performing skeletonization on a movie
%     
% Hunter Elliott
% 12/2010    
%
    properties (SetAccess = protected, GetAccess = public) 
        
        
    end
    
    methods (Access = public)
        
        function obj = SkeletonizationProcess(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = 'Skeletonization';
                super_args{3} = @detectMovieBranches;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                                        
                    funParams.ChannelIndex = 1;                    
                    funParams.MaxRadius = [];%Use the defaults in getBranchesFromMask.m...
                    funParams.SmoothSigma = [];
                    funParams.IsoValue = [];
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'branch_detection'];
                    funParams.RunParallel = 0;
                    funParams.BatchMode = false;                                                      
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
                        
    end
    
end
        
    
    