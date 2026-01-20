function safeExport(target, outFile)
    % target: axes handle or figure handle you're currently drawing into
    if ~isa(target,'matlab.graphics.axis.Axes')
        axSrc = get(target,'CurrentAxes');   % figure given
    else
        axSrc = target;                      % axes given
    end

    % If ancestor is a uifigure, copy content to a classic figure
    isUI = ~isempty(ancestor(axSrc,'matlab.ui.Figure'));
    if isUI
        f = figure('Visible','off','Renderer','opengl','Color','w');
        ax = axes('Parent',f,'Units','normalized','Position',[0 0 1 1]);
        copyobj(allchild(axSrc), ax);                 % copy children only
        % (Optional) sync limits/aspect if needed:
        axis(ax, axis(axSrc)); pbaspect(ax, pbaspect(axSrc));
        set(ax,'XLimMode',axSrc.XLimMode,'YLimMode',axSrc.YLimMode,'CLim',axSrc.CLim);
    else
        f = ancestor(axSrc,'figure'); ax = axSrc;
        set(f,'Renderer','opengl'); set(f,'Color','w');
    end

    exportgraphics(ax, outFile, 'Resolution',300,'BackgroundColor','white');
    if isUI, close(f); end
end