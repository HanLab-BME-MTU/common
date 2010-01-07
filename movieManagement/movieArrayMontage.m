function varargout = movieArrayMontage(movieArray,chanName,iFrame)

% movieArrayMontage(movieArray,chanName,iFrame)
%
% h = movieArrayMontage(movieArray,chanName,iFrame)
%
% Makes a big figure containing an image from each movie in the movieArray
% at the specified channel and frame.
%
% 
% Input:
% 
%   movieArray - A cell-array of movieData structures. The moviedata structures
%   should be formatted as created by setupMovieData.m
% 
%   chanName - A character array containing the name of the channel to
%   display images from.
%   Optional. If not specified, the first channel is displayed.
%
%   iFrame - A positive integer corresponding to the frame number to
%   display images from.
%   Optional. If not specified, the first frame is displayed.
%
%
% Output:
% 
%   h - The handle of the figure the images were displayed on.
% 
% Hunter Elliott
% 
% 3/2009
%

if nargin < 1 || isempty(movieArray)
    error('Come on, you have to input something!')
end

nMovies = length(movieArray(:));

if nargin < 2 || isempty(chanName)
    chanName = movieArray{1}.channelDirectory{1};
end

if nargin < 3 || isempty(iFrame)
    iFrame = 1;
end

%Make the figure
fHan = figure;

%Return the handle if requested
if nargout > 0
    varargout{1} = fHan;
end

%Determine the size of the grid of images
gridSize = ceil(sqrt(nMovies));

%Loop through the movies and display an image from each one.
for j = 1:nMovies
    
    subplot(gridSize,gridSize,j);%Switch to current plot
    axHandle = gca; %Get handle for current axes;
    try
        imageViewer(movieArray{j},'Channel',chanName,'Frame',iFrame,'AxesHandle',axHandle) %Show the image
        title(num2str(j));
    catch
        text(0,0,'Problem with movie channel/frame...','color','w')
    end
end