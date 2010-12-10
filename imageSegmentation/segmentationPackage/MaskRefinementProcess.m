classdef MaskRefinementProcess < MaskProcessingProcess
    %Class definition for post processing using refineMovieMasks.m
    
    
    methods(Access = public)        
        function obj = MaskRefinementProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = 'Mask Refinement';
                super_args{3} = @refineMovieMasks;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%                                            
                    funParams.ChannelIndex = 1:nChan; %Default is to attempt to refine masks for all channels
                    funParams.SegProcessIndex = []; %No default.
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'refined_masks'];      
                    funParams.MaskCleanUp = true;
                    funParams.MinimumSize = 10;
                    funParams.ClosureRadius = 3;
                    funParams.ObjectNumber = 1; %Default is to keep only 1 object per mask
                    funParams.FillHoles = true;
                    funParams.EdgeRefinement = false; %This off by default because it sort of sucks, and is slow.
                    funParams.MaxEdgeAdjust = []; %Use refineMaskEdges.m function defaults for these settings
                    funParams.MaxEdgeGap = [];    
                    funParams.PreEdgeGrow = [];
                    funParams.BatchMode = false;                                              
                end
                %Make sure the input parameters are legit??
                super_args{4} = funParams;                    
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
            
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
            description = 'This process can be used to post-process (and hopefully improve) the masks for selected channels of the movie. This is generally not needed, but can be used in situations where the masks have problems. This is done by removing small background objects, filling any holes in masks, and so on. For more details, see the Settings help for each parameter. The refined masks will be used for the later steps, but will NOT affect the background masks.';
            paramList = {'Input Channels',...
                         'Mask Clean-up',...                     
                         'Minimum Size',...
                         'Closure Radius',...
                         'Object Number',...
                         'Fill Holes',...
                         'Mask Edge Refinement'};                         
                     
            paramDesc = {'This allows you to select which channels you want to refine the masks of. By default, all channels will be refined. Select the channels by clicking on them in the "Available Input Channels" box and then clicking "Select->" to move them to the "Selected Channels" box. You can un-select a channel by clicking the "Delete" button',...
                         'If this box is checked, the basic mask clean-up procedures will be performed. Either this or "Mask Edge Refinement" must be checked. For descriptions of the individual steps in clean-up, see the parameters described below.',...
                         'This parameter specifies the minimum size (in pixels) of objects which should be kept in the mask. Objects smaller than this size will be removed from the mask. This is useful in removing small areas outside the cell which are included in the mask because they are brighter than background, but are not part of a fluorescent object of interest (cell, single molecule, etc).',...
                         'This parameter specifies the size of a disk that will be used to "close" the mask. This has the effect of connecting objects which are closer together than this radius. It also has the effect of rounding sharp corners of the mask whose radius is smaller than this. For additional help, see the help for imclose.m',...
                         'This specifies the number of objects to keep in each mask. That is, if Object Number is set to 2, only the two largest objects in the mask will be kept. This is useful for removing larger background objects, or cells which are only partially in the image, while retaining the cell or object of interest.',...
                         'You guessed it: If this box is checked, any holes in the mask will be filled in. This is useful because sometimes dark areas within a cell (e.g. the nucleus) or other fluorescent object may not be originally included in the mask. Since cells usually don''t have holes in them, filling in these holes is pretty safe.',...
                         'This is very rarely needed, and doesn''t always work. If checked, the edges of the masks will be further refined (This is done after mask clean-up, if it is selected). For more details, and a description of the parameters, see the help section of refineMovieMasks.m'};
                                                  
            if all
                text = makeHelpText(description,paramList,paramDesc);
            else
                text = makeHelpText(description);
            end
             
        end
    end
end
    