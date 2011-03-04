classdef ThresholdProcess < SegmentationProcess
    %A function-specific process for segmenting via thresholding using
    %thresholdMovie.m
    
    methods (Access = public)
        function obj = ThresholdProcess(owner,outputDir, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = 'Thresholding';
                super_args{3} = @thresholdMovie;                           
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'masks'];      
                    funParams.ThresholdValue = []; %Default is to use automatic threshold selection                    
                    funParams.ChannelIndex = 1 : numel(owner.channels_);                                                           
                    funParams.MaxJump = 0; %Default is no jump suppression
                    funParams.GaussFilterSigma = 0; %Default is no filtering.
                    funParams.BatchMode = false;                                              
                end
                %Make sure the input parameters are legit??
                super_args{4} = funParams;                    
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end               
            
    end
end