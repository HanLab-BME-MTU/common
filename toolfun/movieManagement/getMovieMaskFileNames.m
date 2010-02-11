function maskNames = getMovieMaskFileNames(movieData,iChannels)

% maskNames = getMovieMaskFileNames(movieData,iChannel)
%
%This function returns a cell array containing the FULL file name (path AND
%filename) filenames of every mask in each requested channel of the input
%movie. It also verifies that the number of file names returned matches the
%number of images specified in movieData.nImages (number of masks and
%images in each channel must match)
%
% Input:
%
%   movieData - The structure describing the movie, as created by
%   setupMovieData.m
%
%   iChannels - Integer scalar or vector, corresponding to indices of the
%   channel(s) to return mask file names for. Must be input if movie has
%   more than one channel.
%
%
% Output:
% 
%   maskNames - A 1XM cell array, where M is the number of channels,
%   containing an 1xN cell array of the mask file names in that channel,
%   where N is the number of images in that channel.
% 
% Hunter Elliott
% 11/2009

if ~checkMovieMasks(movieData)
    error('Input movieData must have masks specified for at least one channel! Check movieData and masks!')
end

if nargin < 2 || isempty(iChannels)     
    iChannels = find(cellfun(@(x)(~isempty(x) | ischar(x)),movieData.masks.channelDirectory));
    
    if length(iChannels) > 1
        error('Must specify mask channel(s) if more than one channel has masks!')
    end            
end

%Check the specific channels requested
if ~checkMovieMasks(movieData,iChannels)
    error('Problem with masks for specified channel(s)! Check mask files and movieData!!')
end
   
    
nChan = length(iChannels);

maskNames = cell(1,nChan);

for i = 1:nChan
    
    currDir = [movieData.masks.directory filesep movieData.masks.channelDirectory{iChannels(i)}];
    tmpNames = imDir(currDir);    
    
    maskNames{i}  = {tmpNames(:).name};
    
    if length(maskNames{i}) ~= movieData.nImages(iChannels(i))
        error(['Wrong number of masks found in channel "' ...
            movieData.channelDirectory{iChannels(i)} '"' ...
            ' check movieData.nImages and mask directory!'])
    end
    
    maskNames{i} = cellfun(@(x)([currDir filesep x]),maskNames{i},'UniformOutput',false);        
    
end
