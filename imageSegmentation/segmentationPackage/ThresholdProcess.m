classdef ThresholdProcess < SegmentationProcess
    %A function-specific process for segmenting via thresholding using
    %thresholdMovie.m
    
    methods (Access = public)
        function obj = ThresholdProcess(owner,outputDir, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = ThresholdProcess.getName;
                super_args{3} = @thresholdMovie;                           
                
                if nargin < 3 || isempty(funParams)                                       
                    funParams = ThresholdProcess.getDefaultParams(owner,outputDir);
                end
                %Make sure the input parameters are legit??
                super_args{4} = funParams;                    
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end               
            
    end
    methods (Static)
        function name = getName()
            name = 'Thresholding';
        end
        function h = GUI()
            h= @thresholdProcessGUI;
        end
        function methods = getMethods(varargin)
            thresholdingMethods(1).name = 'MinMax';
            thresholdingMethods(1).func = @thresholdFluorescenceImage;
            thresholdingMethods(2).name = 'Otsu';
            thresholdingMethods(2).func = @thresholdOtsu;
            thresholdingMethods(3).name = 'Rosin';
            thresholdingMethods(3).func = @thresholdRosin;            

            ip=inputParser;
            ip.addOptional('index',1:length(thresholdingMethods),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            methods=thresholdingMethods(index);
        end
        
        function funParams = getDefaultParams(owner,outputDir)
            %----Defaults----%
            funParams.OutputDirectory = [outputDir  filesep 'masks'];
            funParams.ThresholdValue = []; % automatic threshold selection
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.MaxJump = 0; %Default is no jump suppression
            funParams.GaussFilterSigma = 0; %Default is no filtering.
            funParams.BatchMode = false;
            funParams.MethodIndx = 1;
        end
    end
        
end