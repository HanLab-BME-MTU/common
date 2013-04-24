function cMap = randomColormap(n,seed)
%RANDOMCOLORMAP produes a random colormap which does not contain grays/white/black
%
% colorMap = randomColormap
% colorMap = randomColormap(n)
% colorMap = randomColormap(n,seed,saturation)
%
% Produces a colormap (an nx3 matrix containing RGB values) for a specified
% number of colors. The colors are of random hue, but always always of
% saturation=1 and value=1 so there are no grays or white or black. By
% specifying an integer seed for the random number generation you can make
% sure you always get the same colormap.
%   
% Input:
%
%   n - number of colors in colormap
%
%   seed - non-negative integer specifying the seed to use for random
%   number generation. Specify to obtain the same colormap each time,
%   otherwise it will be randombly generated each time (the random number
%   generator will be shuffled, but then returned to the state it was in
%   prior to generating the colormap so other processes will not be
%   affected)
%
%   saturation - Amount of saturation to use
%
% Output:
%
% colorMap - an nx3 RGB colormap
%


% Hunter Elliott
% 4/15/2013


if nargin < 1 || isempty(n)
    n = 64;
end

%So we can restore the current state of the random number generator
currState = rng;

if nargin < 2 || isempty(seed)    
    rng('shuffle')
else
    rng(seed)
end


cMap = hsv2rgb([rand(n,1) ones(n,2)]);

%Restore random number generator state
rng(currState);
