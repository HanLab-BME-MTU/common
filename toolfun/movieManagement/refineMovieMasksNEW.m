function movieData = refineMovieMasksNEW(movieData,paramsIn)
% REFINEMOVIEMASKS Performs post-processing to improve masks for an input movie.
% 
% movieData = refineMovieMasks(movieData)
% 
% movieData = refineMovieMasks(movieData,paramsIn)
% 
% 
% This function performs several post-processing steps to refine the masks
% for the input movie, overwriting the existing masks. The available
% post-processing steps are listed below.
% NOTE: THE EXISTING MASKS WILL BE OVERWRITTEN BY THE REFINED MASKS!
% 
% Input:
% 
%   movieData - The structure describing the movie, as created using
%   setupMovieData.m
% 
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('ChannelIndex'-> Positive integer scalar or vector) 
%       The integer indices of the channel(s) to perform background subtraction
%       on. This index corresponds to the channel directories location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels
%
%       ('MaskCleanUp' -> True/False)
%       If true, various operations will be performed to clean up the
%       mask. The operations are listed below, and will be performed in the
%       order they are listed.
%       Optional. Default is True.
%
%       These operations include:
%
%           Clean-up Methods/Parameters:
%
%           ('MinimumSize' -> Positive integer scalar)
%           Any objects in the mask whose area is below this size (in
%           pixels) will be removed via morphological area opening.
%           Optional. Default is 10 pixels. Set to Inf to keep all objects.
%        
%           ('ClosureRadius' -> non-negative integer scalar)
%           If this is parameter greater than zero, the mask will be closed
%           using a disk-shaped structuring element of this radius. This
%           has the effect of connecting previously un-connected components
%           in the mask if they are within 2x this distance of one another.
%           Optional. Default is 3 pixels.
%
%           ('ObjectNumber -> Positive integer scalar)
%           Only this number of the largest objects in the mask will be
%           kept. That is, if this number is 2, only the two largest
%           objects will be kept. This step is performed AFTER the edge
%           refinement (if enabled)
%           Optional. Default is 1. Set to Inf to keep all objects.
%
%           ('FillHoles -> True/False)
%           If true, any holes in any objects in the mask will be filled
%           in.
%           Optional. Default is true.
%
%       ('EdgeRefinement' -> True/False)
%       If true, edge detection will be used to refine the position of the
%       mask edge to the nearest detected edge in the image. This will be
%       done after any of the cleanup procedures listed above.
%       Optional. Default is False.
% 
%           Edge Refinement Parameters:
%
%           NOTE: For descriptions and default values, see refineMaskEdges.m           
%           
%           ('MaxEdgeAdjust' -> Positive Integer scalar)
%           
%           ('MaxEdgeGap' -> Positive integer scalar)
% 
%           ('PreEdgeGrow' -> Positive integer scalar)
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
% 
% Output:
% 
%   movieData - The modified movieData structure, with the mask refinement
%   logged in it, including all parameters used.
% 
%   The refined masks will replace the old masks in each channel's mask
%   directory.
% 
% 
% Hunter Elliott 
% 1/2010
%

%% ----------- Input --------- %%

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end
if nargin < 2
    paramsIn = [];
end

%Make sure the move has been segmented
iSegProc = find(cellfun(@(x)(isa(x,'ThresholdProcess')),movieData.processes_),1);   %TEMPORARY - CHECK FOR SEGPROCESS, NOT JUST THRESHOLDPROCESS HLE
%WHAT IF THERE ARE MULTIPLE SEGMENTATION PROCESSES!!!! TEMP TEMP
if isempty(iSegProc) 
    error('Must create masks before refining masks The input movie does not have any segmentation processes!!')
else
   %Check which channels have foreground masks 
   hasMasks = cellfun(@(x)(~isempty(x)),movieData.processes_{iSegProc}.maskPaths_);
end


%Get the indices of any previous background mask processes from this function                                                                              
iProc = find(cellfun(@(x)(isa(x,'RefineMaskProcess')),movieData.processes_),1);                          

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(RefineMaskProcess(movieData));                                                                                                 
end

nChan = length(movieData.channelPath_);


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);


%----Param Check-----%

if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

if p.ClosureRadius < 0 || ~isequal(round(p.ClosureRadius),p.ClosureRadius)
    error('The closure radius must be a non-negative integer!!')
end

%...TEMPORARY - MORE PARAMETER CHECKS NEEDED!! - HLE

%Make sure all the selected channels have foreground masks.
if any(~hasMasks(p.ChannelIndex))
    warning('Cannot refine masks for some channels, because they do not have masks! \n Please segment these channels before attempting refinement!')
    p.ChannelIndex = p.ChannelIndex(hasMasks);
elseif ~any(hasMasks(p.ChannelIndex))
    error('Cannot refine masks - none of the specified channels have masks!')
end

if ~(p.MaskCleanUp || p.EdgeRefinement)
    error('You must enable mask cleanup and/or edge refinement! Otherwise this function has nothing to do!')
end


%% ----------- Init ----------- %%


maskNames = movieData.processes_{iSegProc}.getMaskFileNames(p.ChannelIndex);

if p.EdgeRefinement %Images are only needed for edge-refinement
    imageNames = movieData.getImageFileNames(p.ChannelIndex);
end

nChan = length(p.ChannelIndex);

if p.ClosureRadius > 0 %If closure is to be performed, create the structuring element
    seClose = strel('disk',p.ClosureRadius(1));    
end

nImages = movieData.nFrames_;   
nImTot = nImages * nChan;


%% ---------------- Mask Refinement --------------- %%


if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, refining masks for channel ' num2str(p.ChannelIndex(1)) ' ...']);     
end        


disp('Starting mask refinement...')

for iChan = 1:nChan

    disp(['Refining masks for channel ' num2str(p.ChannelIndex(iChan)) '...']);
    
    if ~p.BatchMode        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, refining masks for channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end
    
    disp(['Refining masks for for channel # ' num2str(p.ChannelIndex(iChan)) ' : '])    

    currMaskDir = movieData.processes_{iSegProc}.maskPaths_{p.ChannelIndex(iChan)};        
    if p.EdgeRefinement
        currImDir = movieData.channelPath_{p.ChannelIndex(iChan)};
    end
    
    for iImage = 1:nImages;
        
        %Load the mask for this frame/channel
        currMask = imread([currMaskDir filesep maskNames{iChan}{iImage}]);        
        
        % ----- Mask Clean Up ------ %
        
        if p.MaskCleanUp
            
            %Perform initial closure operation
            if p.ClosureRadius > 0
                currMask = imclose(currMask,seClose);            
            end
            
            %Remove objects that are too small
            if ~isinf(p.MinimumSize)                
                currMask = bwareaopen(currMask,p.MinimumSize);                                       
            end

           
        end
        
        
        % --------- Mask Edge-Refinement ------ %
        if p.EdgeRefinement
            
            %Load the current image
            currImage = imread([currImDir filesep imageNames{iChan}{iImage}]);
            
            %Call the edge-refinement function
            currMask = refineMaskEdges(currMask,currImage,...
                p.MaxEdgeAdjust,p.MaxEdgeGap,p.preEdgeGrow);
            
            
        end
        % ---------- Object Selection -------- %
        
        %Keep only the largest objects
        if p.MaskCleanUp && ~isinf(p.ObjectNumber)
                
            %Label all objects in the mask
            labelMask = bwlabel(currMask);

            %Get their area
            obAreas = regionprops(labelMask,'Area');      

            %First, check that there are objects to remove
            if length(obAreas) > p.ObjectNumber 
                obAreas = [obAreas.Area];
                %Sort by area
                [~,iSort] = sort(obAreas,'descend');
                %Keep only the largest requested number
                currMask = false(size(currMask));
                for i = 1:p.ObjectNumber
                    currMask = currMask | labelMask == iSort(i);
                end
            end
        end
        
        % ------ Hole-Filling ----- %
        if p.FillHoles
            currMask = imfill(currMask,'holes');
        end
        
        %Write the refined mask to file, over-writing the previous mask.
        imwrite(currMask,[currMaskDir filesep maskNames{iChan}{iImage}]);
        
        if ~p.BatchMode && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end                 
        
    end
   
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------ Finalize ------ %%


%Store parameters/settings in movieData structure
movieData.processes_{iProc}.setSuccess(true);
movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk

disp('Finished creating background masks!')

