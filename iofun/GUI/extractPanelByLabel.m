function [data, ax, info] = extractPanelByLabel(figFile, labelString)
% DATA: cell array of structs with fields X, Y, Tag, Class
% AX: the axes handle that contained the label
% INFO: struct with axes Position, XLim, YLim for reference

% Open figure invisibly
fig = openfig(figFile, 'invisible');

% Find text object(s) that match the label
tx = findall(fig, 'Type', 'text', 'String', labelString);
if isempty(tx)
    error('No text object with String="%s" found in %s', labelString, figFile);
end

% Get the axes that contain those text objects (some figs have duplicates)
axCandidates = unique(ancestor(tx, 'axes'));
if numel(axCandidates) > 1
    % If more than one, prefer the one in the “top row, middle-ish”.
    % Sort by bottom (desc), then left (asc).
    pos = cell2mat(get(axCandidates, 'Position')); % [left bottom width height]
    [~, order] = sortrows([ -pos(:,2), pos(:,1) ]); % top-most, then left-most
    ax = axCandidates(order(1));
else
    ax = axCandidates;
end

% Collect everything in that axes that looks like XY data
cand = findall(ax);  % all descendants
xyKids = cand(arrayfun(@(h) isprop(h,'XData') && isprop(h,'YData'), cand));

data = {};
for k = 1:numel(xyKids)
    x = get(xyKids(k), 'XData');
    y = get(xyKids(k), 'YData');
    % Flatten cell/NaN-separated segments if needed
    if iscell(x); x = x{1}; end
    if iscell(y); y = y{1}; end
    try
        x = x(:); y = y(:);
        if ~isempty(x) && ~isempty(y) && numel(x)==numel(y)
            s = struct();
            s.X = double(x);
            s.Y = double(y);
            s.Tag = string(get(xyKids(k),'Tag'));
            s.Class = string(class(handle(xyKids(k))));
            data{end+1} = s; %#ok<AGROW>
        end
    catch
        % skip any exotic object that doesn’t behave
    end
end

% Axes info
info = struct('Position', get(ax,'Position'), ...
              'XLim', get(ax,'XLim'), ...
              'YLim', get(ax,'YLim'));

% Optional: close the figure now that we extracted data
close(fig);
end