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
                super_args{2} = 'Bleedthrough Correction';
                super_args{3} = @bleedthroughCorrectMovie;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'bleedthrough_corrected_images'];                      
                    funParams.ChannelIndex = [];%No default
                    funParams.BleedChannelIndex = [];%No default
                    funParams.BleedCoefficients = [];%No default
                    funParams.BatchMode = false;                                                                                
                                     
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

end                                   
            