classdef MSSSegmentationProcess < SegmentationProcess
    %A process for segmenting via multi-scale steerable segmentation
    
    % Sebastien Besson
    
    methods (Access = public)
        function obj = MSSSegmentationProcess(owner,outputDir, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = MSSSegmentationProcess.getName;
                super_args{3} = @getMovieMasksMSS;                           
                
                if nargin < 3 || isempty(funParams)    
                    funParams.ChannelIndex = 1:numel(owner.channels_);
                    funParams.OutputDirectory = [outputDir  filesep 'MSSMasks'];
                    funParams.Scales = [1 2 4]; %Default is no jump suppression
                    funParams.FilterOrder = 3;
                end
                super_args{4} = funParams;                    
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end               
            
    end
    methods
        function sanityCheck(obj)
        end
    end
    methods (Static)
        function name = getName()
            name = 'Multi-scale segmentation';
        end
        function h = GUI()
            h= @mssSegmentationProcessGUI;
        end
    end
        
end