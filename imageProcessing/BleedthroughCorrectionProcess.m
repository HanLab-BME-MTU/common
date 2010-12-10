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
    methods (Static)
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end
            description = '"Bleedthrough" is when the fluorophore from one channel contributes to the intensity in another channel, due to spectral overlap or imperfect filtering. With single-chain (intramolecular) biosensors, the localization of the fluorophores in each channel is always identical, so the bleedthrough contibution is automatically cancelled out during the ratioing step. With dual-chain (intermolecular) biosensors, the localization is almost always at least slightly different for the two fluorophores, meaning that bleedthrough must be corrected. This is accomplished by determining bleedthrough coefficients in a separate experiment where each half of the sensor is imaged independently. These coefficients are then used to correct the actual fluorphore images. The coefficients may be calculated from these bleedthrough experiments using the function calculateMovieBleedthrough.m - see the user''s manual for more information.';
            paramList = {'Input Channel',...
                         'Bleedthrough Channels',...
                         'Bleedthrough Coefficients'};
                         
                         
            paramDesc = {'This allows you to select the channel to bleed-through correct. This should be the activity channel (usually FRET) which will be the numerator in the ratio.',...
                         'This box allows you to select the channels which contain fluorophores which are bleeding through into the channel selected for correction (Usually CFP and YFP)',...
                         'This is where you can input the bleedthrough coefficients that you have determined previously. You must input one coefficient for each bleedthrough channel selected above. Type in the coefficient in the box and then click the "Add->" button. Coefficients can be determined from bleedthrough experiments using the function calculateMovieBleedthrough.m - see the users''s manual for more information.'};
                         
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end    
end                                   
            