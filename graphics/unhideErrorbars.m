function unhideErrorbars(figH)
%UNHIDEERRORBARS shows all hidden errorbars in a figure
%
% SYNOPSIS unhideErrorbars(figH)
%
% INPUT    figH handle of the figure you want to turn the errorbars back
%           on. (complement to hideErrorbars)
%
% c: 05/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find lineHandles in figure
lineHandles = findall(figH,'Type','line');

for lh = lineHandles'
    if findstr(get(lh,'Tag'),'errorBar')
        set(lh,'Visible','on')
    end
end