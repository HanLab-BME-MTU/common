function collectionFigure = collectCurves(figureHandles, recolor, showLegend,varargin)
%COLLECTCURVES will copy curves from several figures into one figure
%
% SYNOPSIS collectCurves(figureHandles, recolor, showLegend, pn/pv)
%
% INPUT    figureHandles : vector of figure handles or cell array with {
%                          figureHandles, pn , pv , ..} , where pn/pv are
%                          property name/ property value pairs to select a
%                          specific subset of graphics objects to collect.
%          recolor       : (opt) [0/{1}/2] whether or not to recolor the curves
%                           0: keep original color
%                           1: each curve has a different color
%                           2: same color for curves from the same figure
%          showLegend    : (opt) [0/{1}] whether or not to show legend
%          pn/pv         : propertyName/propertyValue to be set to each
%                           object, if possible
%
% OUTPUT   collectionFigure : handle to the figure with the collected curves
%
% REMARKS  To better identify curves, it might be helpful to tag them
%           first, e.g. by using plot(x,y,'Tag','myDescription'). After
%           collection, the functions '(un)hideErrorbars' can be very
%           convenient for clarity.
%
%
%c: jonas 04/04
% 11/07 - support for bars (ML7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

%------------
% test input
%------------
if nargin == 0 || isempty(figureHandles) || ((iscell(figureHandles) && ~all(ishandle(figureHandles{1}))) || ( ~iscell(figureHandles) && ~all(ishandle(figureHandles))))
    error('please specify valid figure handles as input for COLLECTCURVES')
end
if nargin < 2 || isempty(recolor)
    recolor = 1;
end
if nargin < 3 || isempty(showLegend)
    showLegend = 1;
end
if nargin < 4 || isempty(varargin)
    varargin = {};
end
%-------------

%-----------------
% make new figure #### - should dock if most are docked
newFH = figure('Name','Collected Curves');
newAxH = axes;

%If figureHandle is a cell array extract handles and pn/pv.
if iscell(figureHandles)
    pnpv = figureHandles(2:end);
    if ~isEven(length(pnpv))
        error('please enter property value for each property name.');
    end
    figureHandles = figureHandles{1};
else
    pnpv = {};
end

% reshape figureHandles so that we can run the loop correctly
figureHandles = figureHandles(:)';

lineCt = 1;
colorCt = 1;
legendCt = 1;

axAndTitle = cell(length(figureHandles),4);
xyz = 'xyz';
%----------------------


%---------------------
% loop through figure handles and copy all lines. Add a tag so that you
% will know where the curves came from
for fh = figureHandles
    
    % to be sure: kill legends
    legH = findall(fh,'Tag','legend');
    if ~isempty(legH)
        %figHadLegend = 1;
        delete(legH);
    else
        %figHadLegend = 0;
    end
    
    % collect figure names and axes labels
    ah = setdiff(findall(fh,'type','axes'),findall(fh,'tag','scribeOverlay'));

    
    fhIdxL = fh==figureHandles;
    for i=1:3
        axAndTitle{fhIdxL,i} = get(get(ah,sprintf('%sLabel',xyz(i))),'String');
    end
    axAndTitle{fhIdxL,4} = get(fh,'Name');
    
    % find lines & reshape. Order backwards, because new line handles are
    % appended to the figure children at the top, and we want to observe
    % the sequence in which the lines were plotted
    lineHandles = findall(ah,'Type','line',pnpv{:});
    % check whether there are bars, and whether any of the lineHandles
    % correspond to the y-axis
    barHandles = findall(ah,'Type','patch',pnpv{:});
    if ~isempty(barHandles)
        % remove the lineHandles whose yData is [0,0]
        badLines = findall(lineHandles,'YData',[0,0]);
        hIdx = ismember(lineHandles,badLines);
        lineHandles(hIdx) = [];
    end
    
    lineHandles = lineHandles(end:-1:1)';
    if ~isempty(lineHandles)
        % loop through the lines, copy to figure and update
        for lh = lineHandles
            
            % copy into new figure
            newH(lineCt) = copyobj(lh,newAxH); %#ok<AGROW>
            
            % change tag
            oldTag = get(newH(lineCt),'Tag');
            if ~isempty(axAndTitle{fh==figureHandles,4})
                newTag = [axAndTitle{fh==figureHandles,4} ' - '  oldTag];
            else
            newTag = ['fig-' num2str(fh) ' ' oldTag];
            end
            set(newH(lineCt),'Tag',newTag);
            
            % change color
            if recolor > 0
                if ~isempty(findstr(oldTag,'errorBar'))
                    % set line color for error bars to same as base line
                    set(newH(lineCt),'Color',extendedColors(colorCt-1+(recolor==2)));
                    
                elseif ~isempty(findstr(oldTag, 'TAfit'))
                    % don't change color
                    
                else
                    % set new color
                    set(newH(lineCt),'Color',extendedColors(colorCt));
                    if recolor == 1
                        colorCt = colorCt+1;
                    end
                end
            end
            
            % try and set pn/pv
            if ~isempty(varargin)
                for pIdx = 1:2:length(varargin)
                    try
                        set(newH,varargin{pIdx},varargin{pIdx+1})
                    end
                end
            end
            
            % collect tags for legend
            if ~isempty(findstr(oldTag,'errorBar')) || ~isempty(findstr(oldTag,'TAfit'))
                % do not add to legend
            else
                legendCell{legendCt,1} = newTag; %#ok<AGROW>
                legendLineH(legendCt,1) = newH(lineCt); %#ok<AGROW>
                legendCt = legendCt + 1;
            end
            
            lineCt = lineCt + 1;
        end
    end
    %
    
    % checked for barHandles above
    barHandles = barHandles(end:-1:1);
    if ~isempty(barHandles)
        for bh = barHandles
            
            % copy into new figure
            newH(lineCt) = copyobj(bh,newAxH); %#ok<AGROW>
            
            % change tag
            oldTag = get(newH(lineCt),'Tag');
            newTag = ['fig-' num2str(fh) ' ' oldTag];
            set(newH(lineCt),'Tag',newTag);
            
            % change color
            if recolor
                if ~isempty(findstr(oldTag,'errorBar'))
                    % set line color
                    set(newH(lineCt),'Color',extendedColors(colorCt-1));
                    
                elseif ~isempty(findstr(oldTag, 'TAfit'))
                    % don't change color
                    
                else
                    % set new color
                    set(newH(lineCt),'FaceColor',extendedColors(colorCt));
                    set(newH(lineCt),'EdgeColor',extendedColors(colorCt));
                    if recolor == 1
                        colorCt = colorCt+1;
                    end
                end
            end
            
            % try and set pn/pv
            if ~isempty(varargin)
                for pIdx = 1:2:length(varargin)
                    try
                        set(newH,varargin{pIdx},varargin{pIdx+1})
                    end
                end
            end
            
            % collect tags for legend
            if ~isempty(findstr(oldTag,'errorBar')) || ~isempty(findstr(oldTag,'TAfit'))
                % do not add to legend
            else
                legendCell{legendCt,1} = newTag; %#ok<AGROW>
                legendLineH(legendCt,1) = newH(lineCt); %#ok<AGROW>
                legendCt = legendCt + 1;
            end
            
            lineCt = lineCt + 1;
        end
    end
    
    % if recoloring per figure, update counter now
    if recolor == 2
        colorCt = colorCt+1;
    end
    
    %     % turn legend back on - does not work for some reason
    %     if figHadLegend
    %         axH = findall(fh,'Type','axes');
    %         legend(axH,'show');
    %     end
end

%---------------------
% now show the legend.
if showLegend
    lh=legend(legendLineH,legendCell);
    lHandles = get(lh);
    if isfield(lHandles,'Interpreter')
        set(lh,'Interpreter','none')
    end
end

%-----------------
% add axes labels, figure name if possible. Use 'all' so that [] becomes
% true
goodLabels = all(cellfun(@(x,y)all(strmatch(x,y)),axAndTitle,repmat(axAndTitle(1,:),size(axAndTitle,1),1)),1);
goodLabels = goodLabels & ~cellfun(@isempty,axAndTitle(1,:));
if goodLabels(1)
    xlabel(newAxH,axAndTitle{1,1});
end
if goodLabels(2)
    ylabel(newAxH,axAndTitle{1,2});
end
if goodLabels(3)
    zlabel(newAxH,axAndTitle{1,3});
end
if goodLabels(4)
    set(newFH,'Name',[get(newFH,'Name'),' - ',axAndTitle{1,4}]);
end

%-----------------------
% nargout if asked for

if nargout > 0
    collectionFigure = newFH;
end