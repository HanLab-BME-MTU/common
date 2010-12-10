function movieData = transformMovieMasks(movieData,paramsIn)
%TRANSFORMMOVIEMASKS Spatially transforms the masks of the input movie
% 
% movieData = transformMovieMasks(movieData)
% 
% movieData = transformMovieMasks(movieData,paramsIn)
% 
% 
% This function performs a spatial transformation on the masks of the
% selected channels of the input movie.
% 
%
% Input:
% 
%   movieData - The movieData object describing the movie, as created using
%   setupMovieDataGUI.m
% 
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       masks to. Masks for different channels will be saved as
%       sub-directories of this directory.
%       If not input, the masks will be saved to the same directory as the
%       movieData, in a sub-directory called "masks"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) to perform mask transformation on. This
%       index corresponds to the channel's location in the array
%       movieData.channels_. If not input, the user will be asked to select
%       from the available channels
%
%       ('SegProcessIndex' -> Positive integer scalar or vector) Optional.
%       This specifies SegmentationProcess(s) to use masks from by its
%       index in the array movieData.processes_; If input as a vector,
%       masks will be used from the process specified by the first element,
%       and if not available for a specific channel, then from the next
%       process etc. If not input, and multiple SegmentationProcesses are
%       present, the user will be asked to select one, unless batch mode is
%       enabled in which case there will be an error.  
%
%       ('TransformFilePaths' -> Cell array of Character strings) A cell
%       array specifying The FULL path and filename of the .mat file
%       containing the transform to apply to the images in each channel.
%       Should contain one element for each channel to be transformed. The
%       transform should be of the format used by imtransform.m. If not
%       input, the user will be asked to locate a file containing a
%       transform for each channel, UNLESS batchmode is enabled, in which
%       case an error will be generated.
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
% 
% Output:
% 
%   movieData - The modified movieData object, with the mask transform
%   logged in it, including all parameters used.
% 
%   The transformed masks will be written to a sub-directory of the
%   OutputDirectory.
% 
% 
% Hunter Elliott 
% 6/2010
%
%% -------- Parameters ---------- %%

dName = 'xformed_masks_for_channel_';%String for naming the mask directories for each channel
pString = 'xformed_mask_'; %Prefix for saving masks to file


%% ----------- Input --------- %%

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end
if nargin < 2
    paramsIn = [];
end


%Get the indices of any previous mask transform processes from this function                                                                              
iProc = movieData.getProcessIndex('MaskTransformationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(MaskTransformationProcess(movieData));                                                                                                 
end


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);
nChanX = numel(p.ChannelIndex);

if isempty(p.SegProcessIndex)    
    if p.BatchMode
        %If batch mode, just get all the seg processes
        p.SegProcessIndex = movieData.getProcessIndex('SegmentationProcess',Inf,0);            
    else
        %Otherwise, ask the user 
        p.SegProcessIndex = movieData.getProcessIndex('SegmentationProcess',1,1);
    end
end

if isempty(p.SegProcessIndex) 
    error('This function requires that the input movie has already been segmented - no valid SegmentationProcesses were found!')
end


nProc = numel(p.SegProcessIndex);
hasMasks = false(nChanX,nProc);
%Check every specified process for masks
for j = 1:nProc

    %Make sure the specified process is a SegmentationProcess
    if ~isa(movieData.processes_{p.SegProcessIndex(j)},'SegmentationProcess')
        error(['The process specified by SegProcessIndex(' num2str(j) ') is not a SegmentationProcess!'])
    end

    %Check which channels have masks from this process
    hasMasks(:,j) = movieData.processes_{p.SegProcessIndex(j)}.checkChannelOutput(p.ChannelIndex);        

end


%Check if transform files have been specified, and if not, get them
for j = 1:nChanX
    
    if isempty(p.TransformFilePaths{p.ChannelIndex(j)})
        [currFile,currPath] = uigetfile('*.mat',...
            ['Please select the transformation file for channel ' ...
            num2str(p.ChannelIndex(j)) ':']);
        
        if currFile == 0
            error('You must specify a transformation file to continue!')
        end        
        p.TransformFilePaths{p.ChannelIndex(j)} = [currPath currFile];                    
    end
    %This method will check validity of file....
    movieData.processes_{iProc}.setTransformFilePath(p.ChannelIndex(j),...
                                                p.TransformFilePaths{p.ChannelIndex(j)});        
end

%Make sure all the selected channels have foreground masks.
if any(~sum(hasMasks,2))
    warning('biosensors:transformMasks:noFGmasks',...
        'Cannot transform masks because some channels do not have foreground masks! Please segment these channels before transforming masks!')
end


%% ----------- Init ----------- %%


%Load the transformations
xForms = movieData.processes_{iProc}.getTransformation(p.ChannelIndex);


%Get original image size. Mask pixels that are transformed out of this
%area will be omitted to preserve this size
m = movieData.imSize_(1);
n = movieData.imSize_(2);

nImages = movieData.nFrames_;   
nImTot = nImages * nChanX;

%Set up the input / output directories
outMaskDirs = cell(1,nChanX);
maskDirs = cell(1,nChanX);
maskNames = cell(1,nChanX);

for j = 1:nChanX;
    
    %Get the first seg process with masks for this channel
    iP = p.SegProcessIndex(find(hasMasks(j,:),1));
    
    movieData.processes_{iProc}.setInMaskPath(p.ChannelIndex(j),...
        movieData.processes_{iP}.outMaskPaths_(p.ChannelIndex(j)));
    
    maskDirs(j) = movieData.processes_{iP}.outMaskPaths_(p.ChannelIndex(j));
    maskNames(j) = movieData.processes_{iP}.getOutMaskFileNames(p.ChannelIndex(j));
    
    %Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))];    
        
    %Check/create directory
    mkClrDir(currDir)               
    
    %Save this in the process object
    movieData.processes_{iProc}.setOutMaskPath(p.ChannelIndex(j),currDir);
    outMaskDirs{j} = currDir;
   
end


%% ---------------- Mask Transformation --------------- %%


if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, transforming masks for channel ' num2str(p.ChannelIndex(1)) ' ...']);     
end        


disp('Starting mask transformation...')

for iChan = 1:nChanX

    disp(['Transforming masks for channel ' num2str(p.ChannelIndex(iChan)) '...']);
    
    if ~p.BatchMode        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, transforming masks for channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end                
            
    for iImage = 1:nImages;
        
        %Load the mask for this frame/channel
        currMask = imread([maskDirs{iChan} filesep maskNames{iChan}{iImage}]);        
        
        %Apply the transformation
        currMask = imtransform(currMask,xForms{iChan},'XData',[1 m],'YData',[1 n],'FillValues',0);
 
        if ~p.BatchMode && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end
        
        %Write the refined mask to file, over-writing the previous mask.
        imwrite(currMask,[outMaskDirs{j} filesep pString maskNames{iChan}{iImage}]);
     
        
    end
   
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------ Finalize ------ %%


%Store parameters/settings in movieData structure

movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk

disp('Finished transforming masks!')

