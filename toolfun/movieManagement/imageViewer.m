function varargout = imageViewer(movieData,varargin)

%
% imageViewer(movieData)
%
% imageViewer(movieData,'OptionName',optionvalue,...)
%
% figHan = imageViewer(movieData)
% 
%
% Displays an image from the input movie data at the specified
% frame number and channel. If no frame number or channel is specified, the
% first image in the first channel (alphabetically) is displayed.
%
% Input: 
% 
%   movieData - the structure describing the movie to view an image from,
%   as created using setupMovieData.m
% 
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('Channel' -> cell array of character array(s))
%       This option should contain a character string representing the name
%       of the channel(s) to display (i.e. the name of the sub-directory for
%       that channel's images)
%       Optional. If not specified, the first channel alphabetically will
%       be displayed.
%
%       ('Frame' -> positive integer)
%       This option specifies the frame number of the image to display.
%       Optional. Default is to display frame 1.
%
%       ('ColorMap' -> matlab colormap)
%       This option specifies a matlab built-in colormap such as jet or
%       cool to use when displaying the image.
%       Optional. Default is gray.
%
%       ('AxesHandle' -> axes handle)
%       The handle of the axes to display the image on.
%       Optional. If not specified, a new figure is created.
%
%       ('Overlay' -> character array)
%       A string containing the name of something to overlay on the image.
%       Optional. If not specified, nothing is overlain.
%
%         Possible overlay options are:
%
%           'Mask' - Overlays the outline of the mask for the displayed
%           image.
%   
%
% Hunter Elliott
% Revamped 1/2010
%

%% ----- Parameters ------ %%

colStr = {'r','g','b'}; %The color to use for overlays.


%% ----------------- Input -------------- %%

if nargin < 1
    error('Must input a movieData structure!')
end

%Check the movieData
movieData = setupMovieData(movieData);

%Parse the optional inputs
[chanName,iFrame,cMap,axHandle,overlayName] = parseInput(varargin);

%----Defaults----%

if isempty(chanName)
    iChan = 1; %Default is to use first channel. This is USUALLY the first alphabetically...
    chanName = movieData.channelDirectory(1);
    nChan = 1;
else    
    if ~iscell(chanName)
        chanName = {chanName};
    end    
    nChan = length(chanName);
    
    if nChan > 3
        error('ImageViewer.m can only display uo to 3 channels simultaneously! Please specify fewer channels!')
    end
    
    iChan = zeros(1,nChan);
    
    for j = 1:nChan    
        %Find the index of the channel(s) specified
        tmp = find(strcmp(chanName{j},movieData.channelDirectory),1);

        if isempty(tmp)
            error(['Movie does not contain a channel named "' chanName{j} ...
                '" ! Check movieData and channel name!'])
        else
            iChan(j) = tmp;
        end
    end
end    

if isempty(iFrame)
    iFrame = 1;%Default is to display first frame
elseif ~isnumeric(iFrame) || round(iFrame) ~= iFrame || any(iFrame > movieData.nImages(iChan))
    error('Invalid frame number specified! Check frame number!')
end

%If no handle given, create figure
if isempty(axHandle)
    figHandle = figure;      
elseif ishandle(axHandle)
    %If axes handle given, get the figure handle for it
    figHandle = get(axHandle,'Parent');
else
    error('Specified axes handle is not a valid handle!')
end



%% -------- Draw the image ---------%%

%Make sure we are on the correct figure & axes
figure(figHandle)
if ~isempty(axHandle)
    set(figHandle,'CurrentAxes',axHandle)
end

%Get file names for images
imNames = getMovieImageFileNames(movieData,iChan);

if nChan > 1
    currImage = zeros([movieData.imSize(:,iChan(1))' 3]);
    for j = 1:nChan
        %load the image    
        currImage(:,:,j) = mat2gray(imread(imNames{j}{iFrame}));
    end
else
    currImage = imread(imNames{1}{iFrame});
end

%Display the image
imshow(currImage,[])

%If the axes didn't exist before, get it's handle now
if isempty(axHandle)
    axHandle = get(figHandle,'CurrentAxes');
end

%If requested, change the colormap
if ~isempty(cMap)
    colormap(cMap)%This has no effect when displaying multiple channels
end



%% ---Draw the overlays---- %%

%Draw the channel name
arrayfun(@(x)(text(20,20*x,chanName{x},'color',colStr{x},'interpreter','none')),1:nChan);

%Draw the frame number
text(20,80,[num2str((iFrame-1)*movieData.timeInterval_s) ' / ' ...
    num2str((movieData.nImages(iChan)-1)*movieData.timeInterval_s) ' s' ],'color',colStr{1},'FontSize',12)

%Get the current axes and turn hold on
hold(axHandle,'on')

if ~isempty(overlayName)
    switch overlayName

        case 'mask'

            hasMasks = arrayfun(@(x)(checkMovieMasks(movieData,x)),iChan);
                
            for j = 1:nChan
                if hasMasks(j)
                    
                    maskNames = getMovieMaskFileNames(movieData,iChan(j));

                    %Load the mask
                    currMask = imread(maskNames{1}{iFrame});

                    %Convert the mask into a boundary
                    maskBounds = bwboundaries(currMask);

                    %Plot the boundar(ies)                
                    cellfun(@(x)(plot(x(:,2),x(:,1),colStr{j})),maskBounds);                
                else                                
                    disp('No valid masks for specified channel(s)! Check movieData.masks and mask directory/files!')                
                end                                    
            end
            
        otherwise
            
            disp(['"' overlayName '" is not a recognized overlay type!'])

    end
end

if nargout > 0
    varargout{1} = figHandle;
end

function [chanName,iFrame,cMap,axHan,overlayName] = parseInput(argArray)

chanName = [];
iFrame = [];
cMap = [];
axHan = [];
overlayName = [];

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
        
        
        case 'Channel'
            chanName = argArray{i+1};
            
        case 'Frame'
            iFrame = argArray{i+1};

        case 'ColorMap'
            cMap = argArray{i+1};
            
        case 'AxesHandle'
            axHan = argArray{i+1};
            
        case 'Overlay'
            overlayName = argArray{i+1};
            
        otherwise

            error(['"' argArray{i} '" is not a valid option name! Please check input!'])
    end
end                       
