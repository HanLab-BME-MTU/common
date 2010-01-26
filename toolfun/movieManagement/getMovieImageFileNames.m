function imNames = getMovieImageFileNames(movieData,iChannels)

% imNames = getMovieImageFileNames(movieData,iChannel)
%
%This function returns a cell array containing the FULL file name (path AND
%filename) filenames of every image in each requested channel of the input
%movie. It also verifies that the number of file names returned matches the
%number of images specified in movieData.nImages
%
% Input:
%
%   movieData - The structure describing the movie, as created by
%   setupMovieData.m
%
%   iChannels - The integer indices of the channel(s) to return image names
%   for. Default is to return names for all channels.
%
%
% Output:
% 
%   imNames - A 1XM cell array, where M is the number of channels,
%   containing an 1xN cell array of the image names in that channel, where
%   N is the number of images in that channel.
% 
% Hunter Elliott
% 11/2009

if ~isfield(movieData,'imageDirectory') || ~isfield(movieData,'channelDirectory')
    error('Invalid movieData input! Check movieData!')
end

if nargin < 2 || isempty(iChannels)       
    iChannels = 1:length(movieData.channelDirectory);
end

nChan = length(iChannels);

imNames = cell(1,nChan);

for i = 1:nChan
    
    currDir = [movieData.imageDirectory filesep movieData.channelDirectory{iChannels(i)}];
    tmpNames = dir([currDir filesep '*.tif']);
    
    %If not .tif found, try .STK
    if isempty(tmpNames)
        tmpNames = dir([currDir filesep '*.STK']);
    end
    
    imNames{i}  = {tmpNames(:).name};
    
    if length(imNames{i}) ~= movieData.nImages(iChannels(i))
        error(['Wrong number of images found in channel "' ...
            movieData.channelDirectory{iChannels(i)} '"' ...
            ' check movieData.nImages and image directory!'])
    end
    
    imNames{i} = cellfun(@(x)([currDir filesep x]),imNames{i},'UniformOutput',false);        
    
end
