function iChannels = selectMovieChannels(movieData,multiple,promptString)

%
% iChannels = selectMovieChannels(movieData,multiple,promptString)
%
% Allows the user to select channel(s) from those available in the input
% movie and returns their index.
%
% Input:
%
%   movieData - Structure describing movie created with setupMovieData.m
%
%   multiple  - If true, multiple channels can be selected, if false, only
%   one can be selected.
%
%   promptString - Character string with the prompt to display above the
%   selection dialogue box.
%
%
% Output:
%
%   iChannels - The integer index of the channels selected.
%
% Hunter Elliott, 10/2009
%

if nargin < 2 || isempty(multiple)
    multiple = true;
end

if nargin < 3 || isempty(promptString)
    if multiple
        promptString = 'Please select the channels:';        
    else
        promptString = 'Please select a channel:';        
    end    
end

if multiple
    sMode = 'multiple';
else
    sMode = 'single';
end


if ~isfield(movieData,'channelDirectory') || length(movieData.channelDirectory) < 1
    error('Invalid movieData! Please check movieData structure!')
end

[iChannels,OK] = listdlg('PromptString',promptString,...
                'SelectionMode',sMode,...
                'ListString',movieData.channelDirectory,...
                'ListSize',[400 500]);

if isempty(iChannels) || ~OK
    error('Must select at least one channel to continue!')
end
        