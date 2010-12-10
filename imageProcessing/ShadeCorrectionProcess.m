classdef ShadeCorrectionProcess < ImageCorrectionProcess
    
    %A class for performing shade correction (illumination non-uniformity
    %correction) on images.
    %
    %Hunter Elliott, 5/2010
    %

    methods (Access = public)
        
        function obj = ShadeCorrectionProcess(owner,outputDir,funParams,shadeImagePaths,...
                                              inImagePaths,outImagePaths)
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = 'Shade Correction';
                super_args{3} = @shadeCorrectMovie;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%      
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'shade_corrected_images'];  
                    funParams.ShadeImageDirectories = []; %No default for this! It will be handled differently...
                    funParams.ChannelIndex = 1:nChan;
                    funParams.MedianFilter = true;
                    funParams.GaussFilterSigma = 0;
                    funParams.Normalize = 1;
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
        
        function sanityCheck(obj)
        % Sanity check will check the correction images
            for i = obj.funParams_.ChannelIndex
                if ~isempty(obj.correctionImagePaths_{i})
                    
                    if ~exist(obj.correctionImagePaths_{i}, 'dir')
                        error('lccb:set:fatal', ...
                            ['The specified shade image channel:\n\n ',obj.correctionImagePaths_{i}, ...
                            '\n\ndoes not exist. Please double check your channel path.'])
                    end
                    fileNames = imDir(obj.correctionImagePaths_{i},true);
                    
                    if isempty(fileNames)
                        error('lccb:set:fatal', ...
                            ['No proper image files are detected in:\n\n ',obj.correctionImagePaths_{i}, ...
                            '\n\nPlease double check your channel path.'])                        
                    end
                                        
                    
                    for j = 1:length(fileNames)
                        imInfo = imfinfo([obj.correctionImagePaths_{i} filesep fileNames(j).name]);
                        if imInfo.Width ~= obj.owner_.imSize_(1) || imInfo.Height ~= obj.owner_.imSize_(2)
                            error('lccb:set:fatal', ...
                                ['Shade correction image - \n\n',...
                                obj.correctionImagePaths_{i},filesep,fileNames(j).name,...
                                '\n\nmust have the same size as input images. Please double check your correction image data.'])
                        end
                    end
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
            description = '"Shade correction" compensates for uneven illumination in the image. The illumination intensity on almost all microscopes higher in the center than near the edges of the movie, often by 15-20%% or more. This illumination pattern also varies with wavelength, laser line, etc. and can therefore introduce error into the final ratio images. This illumination heterogeneity is corrected by taking "shade images" - images with the same illumination and acquisition settings, but taken of a blank area of a coverslip which contains no fluorescent objects. These images are then averaged together and filtered, and used to even-out the illumination in the actual fluorescence images by dividing each image by the illumination pattern.';
            paramList = {'Input Channels',...
                         'Shade Image Channels',...
                         '3x3 Median Filter',...
                         'Gaussian Filter',...
                         'Normalize',...
                         'Normalize to Mean 1',...
                         'Normalize to combined mean'};
                         
            paramDesc = {'This allows you to select which channels you want to perform shade correction on. This should be applied to all channels that are going to be used for ratioing or bleedthrough correction. Select the channels by clicking on them in the "Available Input Channels" box and then clicking "Select->" to move them to the "Selected Channels" box. You can un-select a channel by clicking the "Delete" button',...
                         'This box allows you to specify a directory containing the shade images corresponding to each channel to be corrected. You must specify a directory for each channel to be shade corrected, but the same directory may be specified multiple times (However, this is not recommended!) The directories specified should contain one or more "shade images". It is recommended to take 5 or more, as these will be averaged together to improve the correction. It is also recommended that separate shade-correction images be taken for each channel to be corrected.',...
                         'If this box is checked, a median filter will be applied to the shade images prior to their use as a correction. This is usefull because it minizes the contribution of noise in the shade images, and removes "hot pixels" - pixels which have a much higher-than-normal background value (several hundred counts)',...
                         'If checked, the shade images will also be filtered (smoothed) using a gaussian filter whose sigma (in pixels) is specified by the value in the "Sigma" box. Larger sigmas will give smoother images, but too large of a sigma will cause loss of information in the correction images, and introduce artifacts near the image edge. A good starting value is 1. This can be used to further reduce noise in the shade images and is especially important if only 1 or a small number of shade image(s) are taken.',...
                         'If this box is checked, the shade correction images will be normalized prior to their application. This is highly recommended. Normalization allows the pattern of illumination in the shade images to be used without the actual absolute intensity values affecting the correction.',...
                         'If this option is selected, each shade correction image will be normalized so that its mean is equal to 1 prior to use in correction. This means the images will correct only the spatial pattern of illumination, not the overall illumination intensity.',...
                         'If this option is selected, the shade correction images will be normalized so that the combined mean of all shade images across all channels to be corrected is equal to 1 prior to their use in correction. This allows the relative illumination intensities in the different shade images to be taken into account when correcting. This option is only useful in certain circumstances, and is not generally recommended.'};
            
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end    
end                                   
            