function movieData = makeActivityMovie(movieData,varargin)
%MAKEACTIVITYMOVIE makes and saves a movie of all the images in a specified channel 
% 
% movieData = makeActivityMovie(movieData)
% 
% movieData = makeActivityMovie(movieData,'OptionName',optionValue,...)
% 
% This function creates a .mov or .avi movie of all the images in one
% channel of the input movie.
% 
% Input:
% 
%   movieData - The structure describing the movie, as created with
%   setupMovieData.m
% 
%   ('OptionName'->Possible values)
%
%       ('ChannelIndex'-> positive integer scalar)
%       This number corresponds to the index of the channels to make the
%       movie from.
%       Optional. If not input, user is asked.
%
%       ('FigureHandle' - Positive integer) figure handle to plot the
%       images on.
%       Optional. Default is to create a new figure.
%
%       ('FileName'-> Character array) String specifying file name to save
%       movie as.
%       Optional. Default is "activityMovie_CHANNAME" where CHANNAME will
%       be replaced with the channel's name.
%
%       ('ColorMap' -> name of matlab colormap) The colormap to use to
%       display images for movie making.
%       Optional. Default is jet.
%
%       ('RangeAdjust' -> logical) If true, the colormap range will be
%       adjusted to be the mean +/- 2 std of the first image values,
%       excluding zeros. Otherwise, the colormap range will be span all
%       values in the image.
%       Optional. Default is false.
%
%       ('MakeAvi' - Logical scalar) If true, the movie will be saved as .avi.
%       Optional. Default is false.
%
%       ('MakeMov' - Logical scalar) If true, movie will be saved as .mov.
%       Optional. Default is true.
%
% Output:
%
%   movieData - The movieData with the movie's name and location logged in
%   it.
%
%   The resulting movie will be saved in the analysis directory specified in
%   the movieData.
% 
%Hunter Elliott 
%Re-written 2/2010
%

%% ---------- Input ----------- %%

% Parse input from variable input arguments
[iChannel,figHan,movieName,cMap,rangeAdj,makeAvi,makeMov] = parseInput(varargin);


if isempty(iChannel)
    iChannel = selectMovieChannels(movieData,0,'Select the channel to make a movie of:');
end

if isempty(figHan)
    figHan = fsFigure(.75);
else
    figure(figHan);%Make sure this is the current figure.
end
   
if isempty(movieName)
    movieName = ['activityMovie_' movieData.channelDirectory{iChannel}];
end

if isempty(cMap)
    cMap = 'jet';
end

if isempty(rangeAdj)
    rangeAdj = [];
end

if isempty(makeAvi)
    makeAvi = false;
end

if isempty(makeMov)
    makeMov = true;
end

%% ------ Init ----- %%

imNames = getMovieImageFileNames(movieData,iChannel);
nImages = length(imNames{1});

movieData.movies.activityMovie.status = 0;

%% ------ Movie Making ----- %%


for iImage = 1:nImages
    
    clf(figHan);
    
    currImage = imread(imNames{1}{iImage});
                    
    %Draw the image and configure the axes
    imagesc(currImage)
    colormap(cMap)
    axis image,axis off,axis tight
    colorbar;        
    
    %Adjust the colormap range if requested
    if rangeAdj 
        if iImage == 1                
            %Normalization to first frame for colormap        
            %Find the mean and standard deviation of first image values
            imMean = nanmean(cast(currImage(:),'double'));
            imStd = nanstd(cast(currImage(:),'double'));        
        end
        caxis([imMean-2*imStd imMean+2*imStd])
        
    end
                  
    %Draw the time
    text(10,20,[num2str((iImage-1)*movieData.timeInterval_s) ' s'],'color','w','FontSize',14)

    if makeMov        
        if iImage == 1
            MakeQTMovie('start',[movieData.analysisDirectory filesep movieName '.mov'])
            MakeQTMovie('quality',.7)
        end   
        MakeQTMovie('addfigure')    
    end
    
    if makeAvi
        movieFrames(iImage) = getframe(figHan);
    end
    
end



%% ----- Finalization ---- %%

if makeMov
    MakeQTMovie('finish')
end
if makeAvi
    if isunix
        movie2avi(movieFrames,[movieData.analysisDirectory filesep movieName]);
    else
        movie2avi(movieFrames,[movieData.analysisDirectory filesep movieName],'compression','Cinepak');        
    end
end

movieData.movies.activityMovie.status = 1;
movieData.movies.activityMovie.iFrom = iChannel;
movieData.movies.activityMovie.dateTime = datestr(now);
if makeMov %if both were made, just store the .mov filename
    movieData.movies.activityMovie.fileName = [movieName '.mov'];
else
    movieData.movies.activityMovie.fileName = [movieName '.avi'];
end

if ishandle(figHan)%Make sure the user hasn't closed it already.
    close(figHan);
end


function [iChannel,figHan,mvName,cMap,rangeAdj,makeAvi,makeMov] = parseInput(argArray)
%Sub-function for parsing variable input arguments


%-----Init-----%
iChannel = [];
figHan = [];
mvName = [];
cMap = [];
rangeAdj = [];
makeAvi = [];
makeMov = [];

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
           iChannel = argArray{i+1};
           
       case 'FigureHandle'           
           figHan = argArray{i+1};
   
       case 'FileName'
           mvName = argArray{i+1};
           
       case 'MakeAvi'
           makeAvi = argArray{i+1};
           
       case 'MakeMov'
           makeMov = argArray{i+1};
           
       case 'ColorMap'
           cMap = argArray{i+1};
           
       case 'RangeAdjust'
           rangeAdj = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
               
      
   
end
