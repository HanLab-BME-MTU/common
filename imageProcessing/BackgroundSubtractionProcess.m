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
    methods (Static)
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end
            description = 'Different fluorescence channels and imaging conditions can result in variable background fluorescence - that is, fluorescence in areas where there is no fluorophore. To correct for variations in this background level between channels (and between images) the average background intensity is subtracted from each image. The background area is determined using the background masks, which were created beforehand.';
            paramList = {'Input Channels',...
                         'Background Mask Channels'};                         
                         
            paramDesc = {'This allows you to select which channels you want to perform background subtraction on. This should be applied to all channels that are going to be used for ratioing or bleedthrough correction. Select the channels by clicking on them in the "Available Input Channels" box and then clicking "Select->" to move them to the "Selected Channels" box. You can un-select a channel by clicking the "Delete" button',...
                         'This allows you to select which channels to take background masks from for background subtraction. The default is to use the masks from the same channel that is being correced. However, in the case that one channel has poor quality masks, you may use masks from a different channel to perform background subtraction on that channel. As long as both channels are spatially aligned to within the background mask growth radius, this is a valid alternative.'};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end
end                                   
            