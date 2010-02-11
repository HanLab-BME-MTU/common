function image = imageLoader(movieData,chanID,frameNo)

% 
% image = imageLoader(movieData,chanID,frameNo)
% 
% Loads an image from the movie specified by movieData for the
% specified frame and channel.
% 
% Input:
% 
%     chanID -  Either a string with the channel name, or a positive
%                    Integer specifying the channel number.
%                    Optional. Default is first channel.
%                    
%     frameNo - Frame number of image to load.
%                      Optional. Default is frame 1
%                      
%                      
%  Output:
%  
%     image - the image!!
%     
%     
% Hunter Elliott, 5/2009
% 
% 

if nargin < 2 || isempty(chanID)
    chanID = 1;
end

if nargin < 3 || isempty(frameNo)
    frameNo = 1;
end

movieData = setupMovieData(movieData);

if frameNo > movieData.nImages || round(frameNo) ~= frameNo || frameNo <= 0
    error('Input frame number must be a postive integer less than or equal to the number of images in the specified channel!')
end

if ischar(chanID) %If a channel name was input.
    
    %Check that it matches a channel name
    isChan = strcmpi(chanID,movieData.channelDirectory);
    
    if ~any(isChan)
        error('Input channel name does not match a channel in the input movieData!')
    elseif sum(isChan) > 1
        error('Ambiguous channel names in movieData! Check moviedata!!')
    end
    chanID = find(isChan); %Get the channel index.
    
end

if round(chanID) == chanID && chanID > 0 && chanID <= length(movieData.channelDirectory)
        
    imNames = getMovieImageFileNames(movieData,chanID);
    
    
    if length(imNames{1}) < frameNo
        error('Invalid frame number! Not enough images found!')
    end
    
    image = imread(imNames{1}{frameNo});
    
    
    
else
    error('Input channel ID must be a positive integer specifying the index of a channel in the movieData!')
end
    

    
    
        
        
        
    
