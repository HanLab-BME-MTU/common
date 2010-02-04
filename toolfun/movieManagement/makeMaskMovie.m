function movieData = makeMaskMovie(movieData,varargin)

% movieData = makeMaskMovie(movieData)
% 
% movieData = makeMaskMovie(movieData,'OptionName',OptionValue)
%
% Makes a movie overlaying the outline of the masks on the images from each
% channel so that the user can verify that the segmentation is as desired.
%
% Input:
%
%   movieData - The structure describing the movie created with
%   setupMovieData.m
%
%   Possible Option Names:
%
%   'ChannelIndex' : Positive integer scalar or vector of length <=3
%   These numbers correspond to the indices of the channels to overlay
%   masks from. If not input, user is asked.
%
%   'FigureHandle' - Positive integer figure handle to plot the images / masks on.
%
%   'ShowBackground' - Logical. If true, and background masks have been
%   created using createMovieBackgroundMasks.m, the background mask outline
%   will also be overlain, but as a dotted line.
%
%   'FileName' - Character string specifying file name to save movie as.
%   Default is "maskMovie"
%
%   'MakeAvi' - Logical determining whether the movie is saved as .avi.
%   Default is false.
%
%   'MakeMov' - Logical determining whether or not movie is saved as .mov.
%   Default is true.
%
% Output:
%
%   movieData - The movieData with the movie's name and location logged in
%   it.
%
%   The resulting movie will be saved in the analysis directory specified in
%   the movieData.
% 
% Hunter Elliott, 2/2009
%

%% ------- Input ------------ %%


% Parse input from variable input arguments
[iChannels,figHan,showBkgrnd,mvName,makeAvi,makeMov] = parseInput(varargin);


% ---- Check Input ----%


%Make sure that either makAvi or makeMov are true!!
if ~(makeMov || makeAvi)
    error('Either makeMov or makeAvi MUST be true!')
end

%Check/set up the movieData
movieData = setupMovieData(movieData);

if isempty(figHan)    
    figHan = fsFigure(.75);    
elseif ishandle(figHan)
    figure(figHan)    
else
    error('Input figure handle is not a valid graphics handle!!')
end

if isempty(iChannels)
    iChannels = selectMovieChannels(movieData,true);
end

if ~checkMovieMasks(movieData,iChannels)
    error('Must create masks before making mask movie!!!')    
end

%If the overlay was requested, check the background masks
if showBkgrnd && ~checkMovieBackgroundMasks(movieData,iChannels)
    error('Option showBackground was set to true, but background masks have not yet been created! Create background masks using createMovieBackgroundMasks.m!')
end

%Make sure that all the requested channels have the same number of images
if length(unique(movieData.nImages(iChannels))) > 1
    error('All channels must have the same number of frames to make a mask overlay movie!! Please check the requested channels.')
else
    nImages = movieData.nImages(iChannels(1));
end

nChanMov = length(iChannels);

if nChanMov == 0 || nChanMov > 3
    error('The numer of channels must be between 1 and 3 !!')
end


%% ------ Init ------ %%

%Get the mask directory and file names

maskDir = cell(nChanMov,1);
maskFileNames = cell(nChanMov,1);
for j = 1:nChanMov
    maskDir{j} = [movieData.masks.directory filesep movieData.masks.channelDirectory{iChannels(j)}];
    maskFileNames{j} = dir([maskDir{j} filesep '*.tif']);    
end                       



%Get the image names for each channel
imDir = cell(nChanMov,1);
imNames = cell(nChanMov,1);
for j = 1:nChanMov
    imDir{j} = [movieData.imageDirectory filesep movieData.channelDirectory{iChannels(j)}];
    imNames{j} = dir([imDir{j} filesep '*.tif']);    
end


%If requested, get the background mask directories and names
if showBkgrnd
    bkgrndMaskDir = cell(nChanMov,1);
    bkgrndMaskFileNames = cell(nChanMov,1);
    for j = 1:nChanMov        
        bkgrndMaskDir{j} = [movieData.backgroundMasks.directory filesep ...
                            movieData.backgroundMasks.channelDirectory{iChannels(j)}];
        bkgrndMaskFileNames{j} = dir([bkgrndMaskDir{j} filesep '*.tif']);        
    end
end

%Indicate that the movie-making was started
movieData.masks.movie.status = 0;


%% ------ Movie making -----%%

%Go through the frames/channels and make the movie
mColors = {'r','g','b'};
for j = 1:nImages    
    
    for k = 1:nChanMov
        
        if j == 1 && k == 1
            %Load an image to get the size
            tmp = double(imread([maskDir{k} filesep maskFileNames{k}(j).name]));
            currMasks = zeros([ size(tmp) nChanMov]);
            currImage = zeros([ size(tmp) 3]);            
            if showBkgrnd
                currBackMasks = zeros([ size(tmp) nChanMov]);
            end
        end
                                        
        %Load the mask
        currMasks(:,:,k) = double(imread([maskDir{k} filesep maskFileNames{k}(j).name]));

        %Load the image
        currImage(:,:,k) = mat2gray(double(imread([imDir{k} filesep imNames{k}(j).name])));
        
        if showBkgrnd
            %Load the background masks if requested
            currBackMasks(:,:,k) = double(imread([bkgrndMaskDir{k} filesep bkgrndMaskFileNames{k}(j).name]));
        end
    end
   
    %Show image with mask overlay
    figure(figHan)
    clf
    image(currImage);                
    hold on
    axis off,axis image    
    for k = 1:nChanMov
        maskBorders = bwboundaries(currMasks(:,:,k),4);
        cellfun(@(x)plot(x(:,2),x(:,1),'color',mColors{k}),maskBorders);
        if showBkgrnd
            maskBorders = bwboundaries(currBackMasks(:,:,k),4);
            cellfun(@(x)plot(x(:,2),x(:,1),':','color',mColors{k}),maskBorders);
        end        
    end    
    text(20,20,num2str(j),'color','w','FontSize',16)
    
    
    %Show the channel names and corresponding colors
    for k = 1:nChanMov
        text(max(ylim) - 20,20*k,movieData.channelDirectory{iChannels(k)},'color',mColors{k},'FontSize',16);
    end

            
    if makeAvi
        maskMovie(j) = getframe(figHan); %#ok<AGROW> Pre-allocation not needed with movies. I'm serious.
    end    
           
    if makeMov
        
        if j == 1
            MakeQTMovie('start',[movieData.analysisDirectory filesep mvName '.mov'])
            MakeQTMovie('quality',.85)
        end
        MakeQTMovie('addfigure')
        
    end
end

if makeMov
    MakeQTMovie('finish')
    if ~makeAvi %Only specify the .avi extension if both are made
        movieData.masks.movie.fileName = [mvName '.mov'];
    end
end

if makeAvi
    %Make the movie, save it to the movie directory
    if isunix %Compression is not supported under unix!
        movie2avi(maskMovie,[movieData.analysisDirectory  filesep mvName '.avi']);        
    else
        movie2avi(maskMovie,[movieData.analysisDirectory  filesep mvName '.avi'],'compression','Cinepak')    
    end
    movieData.masks.movie.fileName = [mvName '.avi'];
end


%Modify and save the movieData
movieData.masks.movie.dateTime = datestr(now);
movieData.masks.movie.status = 1;
updateMovieData(movieData);

function [iChannels,figHan,showBkgrnd,mvName,makeAvi,makeMov] = parseInput(argArray)
%Sub-function for parsing variable input arguments



%-----Defaults-----%
iChannels = [];
figHan = [];
showBkgrnd = false;
mvName = 'maskMovie';
makeAvi = false;
makeMov = true;

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName/ value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
              
       case 'ChannelIndex'
           
           iChannels = argArray{i+1};
           
       case 'FigureHandle'
           
           figHan = argArray{i+1};
   
       case 'ShowBackground'
           showBkgrnd = argArray{i+1};
           
       case 'FileName'
           mvName = argArray{i+1};
           
       case 'MakeAvi'
           makeAvi = argArray{i+1};
           
       case 'MakeMov'
           makeMov = argArray{i+1};
           
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
               
      
   
end
