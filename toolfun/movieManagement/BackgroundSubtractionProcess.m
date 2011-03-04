classdef BackgroundSubtractionProcess < ImageCorrectionProcess
    
    %A class for performing background subtraction on images using background masks.
    %
    %Hunter Elliott, 5/2010
    %

    methods (Access = public)
        
        function obj = BackgroundSubtractionProcess(owner,outputDir,funParams,backgroundMaskPaths,...
                                              inImagePaths,outImagePaths)
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = 'Background Subtraction';
                super_args{3} = @backgroundSubtractMovie;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'background_subtracted_images'];                      
                    funParams.ChannelIndex = 1:nChan;                    
                    funParams.MaskChannelIndex = funParams.ChannelIndex;
                    funParams.BatchMode = false;                                                                                

                    
                end
                
                super_args{4} = funParams;    
                
                if nargin > 3           
                    %Set the correction image paths to the background mask paths
                    %input.
                    super_args{7} = backgroundMaskPaths;                
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
        function h = resultDisplay(obj)
           %Overrides default display so subtracted value plots can be
           %shown

           %First, just show the corrected images with the viewer
           h = movieDataVisualizationGUI(obj.owner_,obj);

           %Load and display the averaged correction images
           corrImNames = dir([obj.funParams_.OutputDirectory filesep '*subtraction_values*.mat']);
           if ~isempty(corrImNames)

               for j = 1:numel(corrImNames)
                   figure
                   tmp = load([obj.funParams_.OutputDirectory filesep corrImNames(j).name]);
                   tmpF = fieldnames(tmp);
                   plot(tmp.(tmpF{1}));
                   xlabel('Frame Number')
                   ylabel('Subtracted Background Value, A.U.')
                   title(['Bacground Subtraction Values, Channel ' corrImNames(j).name(end-4) ]);
               end

           end

        end
        
    end

end                                   
            