function figure4ppt(figureH)
%FIGURE4PPT formats a figure for export into powerpoint via copy figure
%
%SYNOPSIS   figure4ppt(figureH)
%
%INPUT      figureH(opt): handle of figure to be formatted (default is current figure)
%
%OUTPUT     the same figure, with:
%               transparent background (gives a crisscrossed pattern on the figure)
%               yellow axes and labels
%               fontsize: 18 for tick marks, 24 for axes labels, 28 for
%               title
%               lineWidth 2
%
%           TO COPY THE FIGURE: goto Edit->Copy Options and select
%           'Preserve Information'. With this, you can even work with the
%           figure as a MicrosoftDrawing Object.
%           CAREFUL: only objects fully visible in the figure are copied
%           Then just copy the figure with Edit->Copy Figure
%
%           The program does not change line colors! However, you can
%           change them easily in Powerpoint with Edit Picture -> Recolor
%
%           DOES NOT WORK PROPERLY UNDER LINUX
%
%c: jonas, 10/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check OS
if ~ispc
    h = warndlg('This does not work properly under linux! (By the way, Powerpoint runs under windows)')
    uiwait(h);
end

%check input
if nargin < 1 | isempty(figureH)
    figureH = gcf;
end

%is it a figure handle?
if ~ishandle(figureH)
    error('input is not a valid handle')
end

%get handles
figureHandles = get(figureH);

%is it really a figure handle
if ~strcmp(figureHandles.Type,'figure')
    error('input is not a figure handle')
end

%set figure background color to none
set(figureH,'Color','none');

%find the axes
childrenH = figureHandles.Children;

%there could be several axes! Therefore loop through all the handles and
%format the axes
for i = 1:length(childrenH)
    
    %do work here only if it is an axis
    if strcmp(get(childrenH(i),'Type'),'axes')
        axesH = (childrenH(i));
        axesHandles = get(axesH);
        
        %background: black, axes: yellow, font: Tahoma18, lineWidth: 2
        set(axesH,'Color','none','FontName','Tahoma','FontSize',18,'LineWidth',2,...
            'XColor', [1 1 0],'YColor', [1 1 0],'ZColor', [1 1 0]);
        
        %Title: tahoma28, rest tahoma24, all yellow
        titleH = axesHandles.Title;
        set(titleH,'FontName','Tahoma','FontSize',28,'Color',[1 1 0]);
        set([axesHandles.YLabel,axesHandles.XLabel,axesHandles.ZLabel],'FontName','Tahoma','FontSize',24,'Color',[1 1 0]);
        
        %all lines: lineWidth 2. No color checks are done!
        lineHList = findall(axesH,'Type','line');
        set(lineHList,'LineWidth',2);
        
    end %if strcmp(get(childrenH(i),'Type'),'axes')
    
end %for i = 1:length(childrenH)

%make figur topmost to show the result
figure(figureH);