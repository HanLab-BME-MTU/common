function varargout = showMovieTimeline(movieData,channelName,nFrames,figureHandle)

% 
% figureHandle = showMovieTimeline(movieData,channel,nFrames,figureHandle)
% 
% This function displays several frames from the input movie side-by-side
% as sub-plots of the same figure.
% 
% Input:
% 
%     channelName - The channel name to view images from, if a
%     multi-channel movie. Default is channel 1.
%         
%     nFrames - The number of frames to display from the movie. Optional,
%     default is 3.
%     
%     figureHandle - The handle of the figure to show the frames on.
%     Optional, if not input, a new figure is created.
%     
% Output:
% 
%     figureHandle - The handle of the figure which the frames were plotted
%     on.
%     
%     
% Hunter Elliott, 10/2009
% 


if nargin < 1 || isempty(movieData)
    error('Must input movieData!')
end

if nargin < 2
    channelName = movieData.channelDirectory{1};
end

if nargin < 3 || isempty(nFrames)
    nFrames = 3;   
end

if nargin < 4 || isempty(figureHandle)    
    figureHandle = figure;
else
    figure(figureHandle)
    clf
    hold on
end

if nFrames <= 3
    nTall = 1;
    nWide = nFrames;    
else
    nTall = ceil(sqrt(nFrames));
    nWide = nTall;
end

iChan = find(strcmp(channelName,movieData.channelDirectory),1);
iFrames = floor(linspace(1,movieData.nImages(iChan),nFrames));

for j = 1:nFrames
    subplot(nTall,nWide,j)
    imageViewer(movieData,'Channel',channelName,'Frame',iFrames(j),'AxesHandle',gca)
end

if nargout > 0
    varargout{1} = figureHandle;
end