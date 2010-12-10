function makeMovieMovie(movieData,varargin)
%MAKEMOVIEMOVIE create .avi or .mov (quicktime) movie of images for the input movie, with optional overlays
%
%   makeMovieMovie(movieData)
%   makeMovieMovie(movieData,'OptionName1',optionValue1,'OptionName2',optionValue2,...)
%
%   This function creates a .avi or .mov (quicktime) file from images in
%   the input movie. The images can be from a single channel or 2-3
%   channels as an RGB overlay. Additionally, certain image analysis
%   results can be overlain on each frame in the output movie. This
%   function uses the function imageViewer.m for image and image overlay
%   display.
% 
%   This function unifies and replaces a bunch of other confusing, shitty
%   functions which served similar purposes such as:
%       makeMaskMovie.m
%       makeActivityMovie.m
%       makeWindowTestingMovie.m
%       makeWindowOverlay.m
%       ....
%
% Input:
% 
%   movieData - The MovieData object describing the movie to make a .avi or
%   .mov from, as created using movieSelectorGUI.m
% 
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option. 
%
%     Possible option name/value pairs:
%
%       ('ConstantScale' -> true/false) If true, the display range (color axis
%       limits) selected for the first image will be used throughout the
%       movie. If false, the range will be selected sparately for each
%       frame.
%       *NOTE*: This option currently only works when a single channel is
%       displayed.
%
%       ('ColorBar' -> true/false) If true, a bar showing the color scale
%       (that is the value associated with each color in the movie) will be
%       displayed. Optional. Default is false (no color bar).
%       NOTE: This option will always be disabled (false) when more than
%       one channel is displayed.
%
%       ('FileName'-> Character array) String specifying file name to save
%       movie as.
%       Optional. Default is "activityMovie"
%
%       ('MakeAvi' -> Logical scalar) If true, the movie will be saved as .avi.
%       Optional. Default is false.
%
%       ('MakeMov' -> Logical scalar) If true, movie will be saved as .mov.
%       Optional. Default is true.             
%
%       **NOTE** Additional options and possible values are all described
%       in the help for imageViewer.m. This includes channel specification,
%       image overlays, colormaps etc.
%
%
% Output:  
%   
%   Each frame of the selected channel(s) will be displayed sequentially,
%   along with any selected overlays. These frames will all be then saved
%   to a .avi and/or .mov file in the movie's specified outputDirectory,
%   with the filename specified by the FileName option
% 
%
%
% Hunter Elliott
% 9/2010
%

%% ------------- Input -------------- %%

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

[constRange,cBar,fName,makeAvi,makeMov,imviewArgs] = parseInput(varargin);


% ---- Defaults ---- %

if isempty(constRange)
    constRange = false;
end
if isempty(cBar)
    cBar = false;
end
if isempty(fName)
    fName = 'activityMovie';
end
if isempty(makeAvi)
    makeAvi = false;
end
if isempty(makeMov)
    makeMov = true;
end


%% ------------- Init ----------------- %%

%Make the figure for display and get the axes handles
figHan = figure;
aHan = gca;

nImages = movieData.nFrames_;

%% ------ Movie Making ----- %%


for iImage = 1:nImages
    
    %Clear the axes so we don't build up a bunch of crap in the figure
    cla(aHan);    
    
    imageViewer(movieData,'Frame',iImage,...
                          'AxesHandle',aHan,imviewArgs{:});
                       
        
    if constRange
        if iImage == 1
            clim = caxis;
        else
            caxis(clim);
        end
    end
    
    if cBar
        %Display a color scale bar
        colorbar
    end
    
    if makeMov        
        if iImage == 1
            MakeQTMovie('start',[movieData.outputDirectory_ filesep fName '.mov'])
            MakeQTMovie('quality',.9)
        end   
        MakeQTMovie('addaxes')    
    end
    
    if makeAvi
        movieFrames(iImage) = getframe(figHan);  %#ok<AGROW> You don't have to initialize with getframe. I promise ;)
    end
    
end



%% ----- Finalization ---- %%

if makeMov
    MakeQTMovie('finish')
end
if makeAvi
    if isunix
        movie2avi(movieFrames,[movieData.outputDirectory_ filesep movieName]);
    else
        movie2avi(movieFrames,[movieData.outputDirectory_ filesep movieName],'compression','Cinepak');        
    end
end


if ishandle(figHan)%Make sure the user hasn't closed it already.
    close(figHan);
end


function [constRange,cBar,fName,makeAvi,makeMov,imviewArgs] = parseInput(argArray)
%Sub-function for parsing variable input arguments



%-----Defaults-----%
constRange = [];
cBar = [];
fName = [];
makeAvi = [];
makeMov = [];
imviewArgs = {};

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
                         
       
       case 'ConstantScale'
           
           constRange = argArray{i+1};
           
       case 'ColorBar'
           cBar = argArray{i+1};
           
       case 'FileName'
           fName = argArray{i+1};
           
       case 'MakeAvi'
           makeAvi = argArray{i+1};
           
       case 'MakeMov'
           makeMov = argArray{i+1};
           
           
       otherwise
                  
           imviewArgs  = [imviewArgs argArray{i:i+1}]; %#ok<AGROW>
           
   end
               
      
   
end
