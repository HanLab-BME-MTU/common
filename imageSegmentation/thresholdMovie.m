function movieData = thresholdMovie(movieData,paramsIn)
%THRESHOLDMOVIE applies automatic or manual thresholding to every frame in input movie
%
% movieData = thresholdMovie(movieData)
% movieData = thresholdMovie(movieData,paramsIn)
%
% Applies manual or automatic thresholding to every frame of the input
% movie and then writes the resulting mask to file as a binary .tif in a
% sub-folder of the movie's analysis directory named "masks"
%
% Input:
% 
%   movieData - A MovieData object describing the movie to be processed, as
%   created by setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
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
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       Optional. The integer index of the channel(s) to segment. If not
%       input, all channels will be segmented. If multiple
%       channels are selected, masks are generated for each independently.
%
%       ('ThresholdValue' -> Positive integer scalar or vector) Optional.
%       Intensity value to threshold images at. If not input, a value is
%       automatically calculated based on the image intensity histogram. If
%       a scalar, the same value is used for all channels. If a vector,
%       each element specifies the value to use for a specific channel.
%       Must be the same order and size as ChannelIndex. (requires good
%       SNR, significant area of background in image).
% 
%       ('MaxJump' -> Positive scalar) If this is non-zero, any changes in
%       the auto-selected threshold value greater than this will be
%       suppressed by using the most recent good threshold. That is, if
%       MaxJump is set to 2.0 and the threshold changes by a factor of 2.2
%       between two consecutive frames, the new threshold will be ignored
%       and the previous threshold used. This option is ignored if the user
%       specifies a threshold.
%       Optional. Default is 0 (no jump suppression)
%
%       ('GaussFilterSigma' -> Positive scalar) If this is entered, the
%       image will be filtered first before thresholding and segmentation.
%       The filter kernel will be taken as a Gaussian with the input sigma.
%       Optional. Default is 0 (no filtering)
%           -- KJ
% 
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - the updated MovieData object with the thresholding
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The masks are written to the directory specified by the parameter
%   OuptuDirectory, with each channel in a separate sub-directory. They
%   will be stored as binary, bit-packed, .tif files. 
%
%
% Hunter Elliott, 11/2009
% Revamped 5/2010
%
%% ----- Parameters ----- %%

pString = 'mask_'; %Prefix for saving masks to file
pfName = 'threshold_values_for_channel_'; %Prefix for saving threshold values to file. Actual file name will have channel number appended.
dName = 'masks_for_channel_';%String for naming the mask directories for each channel

%% ----------- Input ----------- %%


%Check that input object is a valid moviedata TEMP
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end


%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('ThresholdProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ThresholdProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

nChan = numel(movieData.channels_);

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end


nChanThresh = length(p.ChannelIndex);
if ~isempty(p.ThresholdValue)
    if length(p.ThresholdValue) == 1
        p.ThresholdValue = repmat(p.ThresholdValue,[1 nChanThresh]);
    elseif length(p.ThresholdValue) ~= nChanThresh
        error('If you specify a threshold value, you must either specify one value to use for all channels, or 1 value per channel!');    
    end
end



%% --------------- Init ---------------%%

disp('Starting thresholding...')

%Set up the mask directories as sub-directories of the output directory
for j = 1:nChanThresh;
    
    %Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))];    
    %Save this in the process object
    movieData.processes_{iProc}.setOutMaskPath(p.ChannelIndex(j),currDir);
   
    %Check/create directory
    mkClrDir(currDir)               
end

allThresholds = cell(nChanThresh,1);

imageFileNames = movieData.getImageFileNames(p.ChannelIndex);


nImages = movieData.nFrames_;   
nImTot = nImages * nChanThresh;

%Get mask and image directories
maskDirs  = movieData.processes_{iProc}.outMaskPaths_(p.ChannelIndex);
imDirs  = movieData.getChannelPaths(p.ChannelIndex);
    

%% ----- Thresholding ----- %%

if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, thresholding channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end        


for iChan = 1:nChanThresh
        
        
    if ~p.BatchMode        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, thresholding channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end        
    disp(['Thresholding images for channel # ' num2str(p.ChannelIndex(iChan)) ' : '])
    disp(imDirs{iChan})
    disp('Masks will be stored in directory :')
    disp(maskDirs{iChan})
    
    %Initialize vector for threshold values
    allThresholds{iChan} = nan(nImages,1);
    
    for iImage = 1:nImages
        

        %Load the current image        
        currImage = imread([imDirs{iChan} filesep imageFileNames{iChan}{iImage}]);

        %KJ: filter image before thesholding if requested
        if p.GaussFilterSigma > 0
            currImage = Gauss2D(double(currImage),p.GaussFilterSigma,1);
        end

        
        if isempty(p.ThresholdValue)
            try
                currThresh = thresholdFluorescenceImage(currImage);             
            catch %#ok<CTCH>
                %If auto-threshold selection fails, and jump-correction is
                %enabled, force use of previous threshold
                if p.MaxJump > 0
                    currThresh = Inf;
                else
                    if ~p.BatchMode && ishandle(wtBar)
                        warndlg(['Could not automatically select a threshold in frame ' ...
                        num2str(iImage) '! Try specifying a threshold level, or enabling the MaxJump option!']);
                        close(wtBar)
                    end                                                
                    error(['Could not automatically select a threshold in frame ' ...
                        num2str(iImage) '! Try specifying a threshold level, or enabling the MaxJump option!']);
                    
                        
                end
            end
        else            
            currThresh = p.ThresholdValue(iChan);
        end
        
        if p.MaxJump > 0
            %Check the threshold
            if iImage == 1
                allThresholds{iChan}(iImage) = currThresh; %Nothing to compare 1st frame to
            else
                if abs(currThresh / allThresholds{iChan}(find(~isnan(allThresholds{iChan}),1,'last'))-1) > p.MaxJump
                    %If the change was too large, don't store this threshold
                    %and use the most recent good value
                    allThresholds{iChan}(iImage) = NaN;
                    currThresh = allThresholds{iChan}(find(~isnan(allThresholds{iChan}),1,'last'));
                else
                    allThresholds{iChan}(iImage) = currThresh;
                end 
            end
        else
            allThresholds{iChan}(iImage) = currThresh;
        end
        
        %Apply the threshold to create the mask
        imageMask = currImage > currThresh;
    
        %write the mask to file                    
        imwrite(imageMask,[maskDirs{iChan} filesep pString imageFileNames{iChan}{iImage}]);
        
        if ~p.BatchMode && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end
                
    
    end    
   
    
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------ Finish - Save parameters and movieData ----- %%


%Save the threshold values to the analysis directory as seperate files for
%each channel
for i = 1:nChanThresh    
    thresholdValues = allThresholds{i}; %#ok<NASGU>
    save([p.OutputDirectory filesep pfName num2str(p.ChannelIndex(i)) '.mat'],'thresholdValues');
end



movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk


disp('Finished thresholding!')




