function [out] = multipleProgressText(text_, nStep_)
%multipleProgressText shows progress of a loop as text on the screen. Can handle multiple levels of progress
%
%SYNOPSIS
%   [] = multipleProgressText(text_, nStep_)
%       Initializes the progressText
%   [] = multipleProgressText(text_)
%       Changes the optional text and updates progressText
%   [] = multipleProgressText('reset')
%       Resets the program
%   [] = multipleProgressText()
%       Updates progressText
%
%INPUT
%   text        : Optional text in the progress display
%   nSteps      : Max number of steps until completion
%
%OUTPUT
%   out     : diagnostic information about the persistent variables
%       .level      : The layer / level of prgoress display. Larger the
%                     level, smaller the increase in progress fraction
%       .iStep      : array of progress on each level
%       .nStep      : array of max number of steps needed to complete each
%                     level
%Tae H Kim, July 2015

%% Initialization
%persistent
persistent text iStep nStep frac weight level
if isempty(level)
    level = 0;
end

%% Input
if nargin == 0
    iStep(level) = iStep(level) + 1;
end
if nargin == 1 && level == 1
    iStep(level) = iStep(level) + 1;
    text = text_;
end
if nargin == 1 && strcmp(text_, 'reset')
    level = 0;
    iStep = [];
    nStep = [];
    frac = [];
end
if nargin == 2
    level = level + 1;
    iStep(level) = 0;
    nStep(level) = nStep_;
    weight(level, 1) = 1 / prod(nStep(2:end));
    if level == 1
        text = text_;
    end
end
%input check
%throws error if level < 1;
if level < 1
    error('Progress display error: number of steps do not match the number times the program was called -or- program was never initialized');
end

%% Fraction calculation
frac(level) = iStep(level) / nStep(level);

%% Progress Display
if iStep(1) < nStep(1) && iStep(level) ~= nStep(level) && nargin ~= 2
    progressText(frac * weight, text);
elseif iStep(1) == nStep(1)
    progressText(1, text);
end

%% level check
if iStep(level) == nStep(level)
    level = level - 1;
    frac = frac(1:end-1);
    weight = weight(1:end-1);
end

%% Output
if nargout > 0
    out.level = level;
    out.iStep = iStep;
    out.nStep = nStep;
end

end

