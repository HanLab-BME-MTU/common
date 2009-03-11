function progressText(fractionDone,text)
%PROGRESSTEXT shows progress of a loop as text on the screen
%
% SYNOPSIS: progressText(fractionDone,text)
%
% INPUT fractionDone: fraction of loop done (0-1)
%		text (opt): {yourTexthere} : XX% done xx:xx:xx remaining
%                   Note: text can be changed for every call
% OUTPUT none
%
% EXAMPLE
%   n = 1000;
%   progressText(0,'Test run') % Create text
%   for i = 1:n
%       pause(0.01) % Do something important
%       progressText(i/n) % Update text
%   end
%
% REMARKS progressText will set lastwarn to ''
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn based on progressbar.m
% DATE: 29-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent starttime lastupdate clearText printText finalText warnText

% constants
nCharsBase = 27; % change this if changing output format

% Test input
if nargin < 1 || isempty(fractionDone) 
    fractionDone = 0;
end
if nargin < 2 || isempty(text)
    text = '';
else
    text = [text,' : '];
end

if fractionDone == 0 || isempty(starttime)
    % set up everything
    
    % Set time of last update to ensure calculation
    lastupdate = clock - 1;
    
    % Task starting time reference
    starttime = clock;
    
    % create fprintf-expression
    
    printText = sprintf('%s%%2d%%%% done %%s remaining',text);
    initialText = sprintf('%s 0%%%% done xx:xx:xx remaining',text);
    finalText = sprintf('%s100%%%% done %%s elapsed\n',text);
    % get length of fprintf expression
    nChars = nCharsBase + length(text);
    clearText = repmat('\b',1,nChars);
    
    % print initialText and return
    fprintf(1,initialText);
    
%     %fprintfExpression removes old expression before overwriting
%     fprintfExpression = [clearText printText];
%     fprintfExpressionFinal = [clearText, finalText];

% empty warning
lastwarn('');
warnText = '';
    
    return
elseif ~isempty(text)
    % text has been changed. Create fprintfExpressions first, then update
    % clearText
    printText = sprintf('%s%%2d%%%% done %%s remaining',text);
    finalText = sprintf('%s100%%%% done %%s elapsed\n',text);
    fprintfExpression = [clearText printText, warnText];
    fprintfExpressionFinal = [clearText, finalText, warnText];
    
    nChars = nCharsBase + length(text);
    clearText = repmat('\b',1,nChars);
% elseif ~isempty(lastwarn)
%     % add warnings to the end of the progressText
%     % find warning
%     w = lastwarn;
%     nw = length(w);
%     % erase warning
%     fprintf(1,repmat('\b',1,11+nw));
%     % create new warnText
%     w = regexprep(w,sprintf('(%s\s+)',char(10)),' - ');
%     warnText = [warnText,sprintf('\n%%3d - Warning : %s',)
else
    % all is normal. Just generate output
    fprintfExpression = [clearText printText, warnText];
    fprintfExpressionFinal = [clearText, finalText, warnText];
end

% write progress
percentDone = floor(100*fractionDone);

% get elapsed time
runTime = etime(clock,starttime);

% check whether there has been a warning since last time
if ~isempty(lastwarn)
    lastwarn('');
    fprintfExpression = regexprep(fprintfExpression,'(\\b)*','\\n');
    fprintfExpressionFinal = regexprep(fprintfExpressionFinal,'(\\b)*','\\b\\n');
end

if percentDone == 100 % Task completed
    fprintf(1,fprintfExpressionFinal,convertTime(runTime)); % finish up
    clear starttime lastupdate clearText printText finalText % Clear persistent vars
    return
end

% only update if significant time has passed
if etime(clock,lastupdate) < 0.3
    return
end

% update
timeLeft = runTime/fractionDone - runTime;
fprintf(1,fprintfExpression,percentDone,convertTime(timeLeft));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfcn
function timeString = convertTime(time)

timeStruct = sec2struct(time);
if timeStruct.hour > 99
    timeString = '99:59:59';
else
    timeString = timeStruct.str(1:end-4);
end
