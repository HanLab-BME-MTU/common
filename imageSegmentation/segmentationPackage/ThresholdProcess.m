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
    methods (Static)
        function text = getHelp(all)
            %Note: This help is designed for the GUI, and is a simplified
            %and shortened version of the help which can be found in the
            %function.
            if nargin < 1  % Static method does not have object as input
                all = false;
            end
            description = 'This process creates masks for the movie which seperate objects (e.g. fluorescent cells, single fluorophores etc.) from the background. Masks are binary images which contain 1 where there is an object of interest (cell), and zero where there is background. Thresholding creates these masks by finding areas of the image which are brighter than a specific value, called the threshold. These will be created for any channels the user selects, or for all channels by default. The masks can be checked by clicking the "Result" button, which will show the outline of the mask on top of the fluorescence image.';
            paramList = {'Input Channels',...
                         'Use Automatic Thresholding',...
                         'Fixed Threshold',...
                         'Maximum Threshold Jump'};
                         
            paramDesc = {'This allows you to select which channels you want to create masks for. Generally, this is only done for fluorescence channels, and DIC/Phase images are not thresholded. Select the channels by clicking on them in the "Available Input Channels" box and then clicking "Select->" to move them to the "Selected Channels" box. You can un-select a channel by clicking the "Delete" button',...
                         'Select this box to allow the software to attempt to automatically select a threshold level for each image in each selected channel. Depending on the data, this automatic thresholding may fail in some cases. In this case, manual thresholding can be used.',... 
                         'If you un-check "Use automatic thresholding", then this option allows you to manually specify an an intensity value to use for thresholding for each channel. This value should be just barely above the average backround intensity in the image. You can use the matlab function imtool, or ImageJ to determine a threshold to use.',...
                         'This option should be used if the automatic thresholding works on most frames of the movie but fails on a few. If selected, this sets an upper limit on the change of the threshold between frames. For example, if this was set to .5 and the automatically-selected threshold changed by a factor of .51 between two consecutive frames, this new threshold would be ignored and the last good threshold value used.'};
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end
end