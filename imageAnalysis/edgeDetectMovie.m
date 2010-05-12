function movieData = edgeDetectMovie(movieData,varargin)
%EDGEDETECTMOVIE uses matlab built in edge detection on all frames of movie 
% 
% movieData = edgeDetectMovie(movieData)
% movieData = edgeDetectMovie(movieData,'OptionName',optionValue)
% 
% This function simply applies the matlab edge-detection function edge.m
% (using the Canny method) to detect edges in every frame of the input
% movie. The detected edges are saved as binary masks. To convert these to
% masks covering an entire object which contains edges, use
% refineMovieMasks.m after this function.
% 
% 
% Input:
% 
%   movieData - The structure describing the movie, as created using
%   setupMovieData.m
% 
% 
%   Possible Options:
%       ('OptionName' -> possible values)
% 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) 
%       The integer indices of the channel(s) to perform edge-detection
%       on. This index corresponds to the channel directories location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels
% 
%       ('FilterRadius'->Positive scalar)
%       The radius, in pixels, of the filter used to detect edges in the
%       images.
%       
%       ('ThresholdValues'->Positive 1x2 vector)
%       This option specifies the threshold values to use for Canny edge
%       detection. The first value is the low threshold and the secon is
%       the high threshold.
%       Optional. If not specified, automatic threshold selection will be
%       used.
%
%       ('ThresholdAdjust' -> Postive scalar / 1x2 Vector)
%       If specified, the initial threshold which edge.m automatically
%       selects will be multiplied by this number, and a second round of
%       edge-detection will be performed. For Canny edge detection, there
%       are two thresholds and if ThresholdAdjust is a vector, the low
%       threshold will be multiplied byt ThreholdAdjust(1) and the High
%       threshold by ThresholdAdjust(2)
%       Note: This option only applies if a ThresholdValues has NOT been
%       specified.
% 
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
%
% Output:
% 
%   The detected edges are saved as masks to the movie's mask directory as
%   binary .tif files.
%
%
% Hunter Elliott
% 4/2010
%

%% ---- Parameters --- %%

pString = 'mask_'; %Prefix for saving masks to file


%% ------------ Input ----------- %%

movieData = setupMovieData(movieData,'masks');

[batchMode,iChannels,filtRad,threshVals,threshAdj] = parseInput(varargin);

if isempty(batchMode)
    batchMode = false;
end

if isempty(iChannels)    
    if ~batchMode
        iChannels = selectMovieChannels(movieData,1,'Select the channels to edge-detect:');                
    else        
        error('In batch mode, you must specify the option ChannelIndex!')
    end
end

%Determine number of channels
nChan = length(iChannels);

if isempty(filtRad)
    filtRad = 1;
end

%If a single adjust value was given, use it for both thresholds
if ~isempty(threshAdj) && length(threshAdj) == 1
    threshAdj = ones(1,2) * threshAdj;
end


   
%% -------- Edge Detection ---------- %%



imageFileNames = getMovieImageFileNames(movieData,iChannels);

for iChan = 1:nChan
          
    %Get/create directory for masks
    maskDir = [movieData.masks.directory filesep movieData.channelDirectory{iChannels(iChan)}];    
    if ~exist(maskDir,'dir')
        mkdir(maskDir);
    end             
    
    if ~batchMode
        wtBar = waitbar(0,['Please wait, detecting edges for channel "' ...
            movieData.channelDirectory{iChannels(iChan)} '"...']);
    end
            
    for iImage = 1:movieData.nImages(iChannels(iChan));
            
        currIm = imread(imageFileNames{iChan}{iImage});
        
        if isempty(threshVals) && ~isempty(threshAdj)
            [~,currThresh] = edge(currIm,'canny',[],filtRad);
            currThresh = currThresh .* threshAdj;
        elseif ~isempty(threshVals)            
            currThresh = threshVals;
        else
            currThresh = [];
        end
        
        currEdges = edge(currIm,'canny',currThresh,filtRad);
        
        %write the mask to file
        currFname = imageFileNames{iChan}{iImage}...
            (max(regexp(imageFileNames{iChan}{iImage},filesep))+1:end);
        imwrite(currEdges,[maskDir filesep pString currFname]);
        
        if ~batchMode
            waitbar(iImage/movieData.nImages(iChannels(iChan)),wtBar)
        end                
    
    end
                  
end

if ~batchMode && ishandle(wtBar)
    close(wtBar)
end
    
%% ------ Finalize ----- %%


movieData.masks.status = 1;
movieData.masks.dateTime = datestr(now);
movieData.masks.parameters.threshAdj = threshAdj;
movieData.masks.parameters.filterRadius = filtRad;
movieData.masks.iFrom(iChannels) = iChannels; %The from and to channel indices are the same.
movieData.masks.channelDirectory(iChannels) = movieData.channelDirectory(iChannels);

updateMovieData(movieData);

disp('Finished edge detection!')
    

function [batchMode,iChannels,filtRad,threshVals,threshAdj] = parseInput(argArray)


%Init output

batchMode = [];
iChannels = [];
filtRad = [];
threshVals = [];
threshAdj = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
              
       case 'BatchMode'
           batchMode = argArray{i+1};
           
       case 'ChannelIndex'
           iChannels = argArray{i+1};
           
       case 'FilterRadius'
           filtRad = argArray{i+1};
           
       case 'ThresholdAdjust'
           threshAdj = argArray{i+1};
           
       case 'ThresholdValues'           
           threshVals = argArray{i+1};
           
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end