classdef BleedthroughCorrectionProcess < ImageCorrectionProcess
    
    %A class for performing bleedthrough correction on images.
    %
    %Hunter Elliott, 5/2010
    %

    methods (Access = public)
        
        function obj = BleedthroughCorrectionProcess(owner,outputDir,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                
                
                super_args{1} = owner;
                super_args{2} = BleedthroughCorrectionProcess.getName;
                super_args{3} = @bleedthroughCorrectMovie;                               
                if nargin<2, outputDir = owner.outputDirectory_; end
                if nargin < 3 || isempty(funParams) 
                    
                    funParams=BleedthroughCorrectionProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;    
                
                if nargin > 3           
                    %Set the correction image paths to the shadeImage paths
                    %input.
                    super_args{7} = shadeImagePaths;                
                end
                
                if nargin > 4
                    super_args{5} = inImagePaths;
                end
                
                if nargin > 5
                    super_args{6} = outImagePaths;
                end                
                
            end
            
            obj = obj@ImageCorrectionProcess(super_args{:});
        end   
        
    end
    methods (Static)
        function name =getName()
            name = 'Bleedthrough Correction';
        end
        function h = GUI()
            h= @bleedthroughCorrectionProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'bleedthrough_corrected_images'];
            funParams.ChannelIndex = [];%No default
            funParams.ProcessIndex = [];%No default
            funParams.BleedChannelIndex = [];%No default
            funParams.BleedCoefficients = [];%No default
            funParams.BatchMode = false;      
        end
    end
end                                           