classdef TracksDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Linestyle='-';
        Linewidth=1;
        GapLinestyle='--';
        Color = 'r';
        MergeColor = 'y';
        MergeMarker = 's';
        SplitColor = 'g';
        SplitMarker = 's';
        useDragtail=true;
        dragtailLength=10;
        showLabel=false;
        markMergeSplit=false;
        ButtonDownFcn=[];
    end
    methods
        function obj=TracksDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj, tracks, tag, varargin)
            if isempty(tracks), h = -1; return; end
            % Get track length and filter valid tracks
            trackLengths = cellfun('prodofsize',{tracks.xCoord});
            validTracks = trackLengths>0;
            tracks = tracks(validTracks);
            trackLengths = trackLengths(validTracks);
            
            nTracks = numel(tracks);
            
            % Constraing the dragtail length between 2 and the maximum
            % track length
            if obj.useDragtail
                dLength = max(2,min(obj.dragtailLength,max(trackLengths)));
            else
                dLength = max(trackLengths);
            end
            
            % Concatenate data in a matrix of size dragtailLength x nTracks
            xData = NaN(dLength, nTracks);
            yData = NaN(dLength, nTracks);
            for i = 1 : max(trackLengths)
                selected = trackLengths == i;
                xTemp = vertcat(tracks(selected).xCoord)';
                yTemp = vertcat(tracks(selected).yCoord)';
                if(i < dLength)
                    xData(1:i,selected) = xTemp;
                    yData(1:i,selected) = yTemp;
                else
                    xData(:,selected) = xTemp(end-dLength+1:end,:);
                    yData(:,selected) = yTemp(end-dLength+1:end,:);
                end
            end
            displayLength = trackLengths;
            displayLength(trackLengths > dLength) = dLength;
            
            % Initialize matrix for gaps
            xGapData = NaN(size(xData));
            yGapData = NaN(size(xData));
            
            % Label gaps: series of NaNs not-connected to the border
            I = isnan(xData);
            I = [I; zeros(size(I))];
            I = reshape(I, size(I,1)/2, size(I,2)*2);
            I = [zeros(size(I,1), 1) I];
            I = imclearborder(I);
            
            cc = bwconncomp(I);
            gapLengths = cellfun('length',cc.PixelIdxList);

            % Obtain first and last coordinate indices
            iFirstLast = zeros(2,cc.NumObjects);
            uGapLengths = unique(gapLengths);
            for i=uGapLengths(:)'
                s = i == gapLengths;
                temp = [cc.PixelIdxList{s}];
                iFirstLast(:,s) = temp([1 i],:);
            end

            % Remove extra columns by adjusting index
            m = mod(iFirstLast,cc.ImageSize(1));
            iFirstLast = (iFirstLast + m - cc.ImageSize(1))/2;
            
            
            iFirst = iFirstLast(1,:) - 1;
            iLast = iFirstLast(2,:) + 1;
            xFirst = xData(iFirst);
            xLast = xData(iLast);
            yFirst = yData(iFirst);
            yLast = yData(iLast);
            
            % Fill gaps
            for gapLength = 1 : max(gapLengths)
                s = gapLengths == gapLength;
                px = [cc.PixelIdxList{s}];
                
                % Remove extra columns by adjusting index
                m = mod(px,cc.ImageSize(1));
                px = (px + m - cc.ImageSize(1))/2;
                
                a = (1:gapLength)'/(gapLength+1);
                xGapData(iFirst) = xFirst;
                xGapData(iLast) = xLast;
                yGapData(iFirst) = yFirst;
                yGapData(iLast) = yLast;
                xGapData(px) = ...
                    bsxfun(@plus,bsxfun(@times,a,xLast(s) - xFirst(s)),xFirst(s));
                yGapData(px) = ...
                    bsxfun(@plus,bsxfun(@times,a,yLast(s) - yFirst(s)),yFirst(s));
            end
            
            eventsExist = isfield(tracks,'splitEvents') || isfield(tracks,'mergeEvents');

            % Initialize matrix for split events
            if(eventsExist)
                hasSplitEvents = ~cellfun('isempty',{tracks.splitEvents})';
            end
            xSplitData = NaN(dLength, nTracks);
            ySplitData = NaN(dLength, nTracks);
            for i = find(hasSplitEvents)'
                eventTimes = false;
                eventTimes([tracks(i).splitEvents tracks(i).splitEvents+1]) = true;
                eventTimes = find(eventTimes);
                dragtailWindow = [trackLengths(i) - displayLength(i) + 1 trackLengths(i)];
                eventTimes = eventTimes(eventTimes >= dragtailWindow(1) & eventTimes <= dragtailWindow(2));
                xSplitData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).xCoord(eventTimes);
                ySplitData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).yCoord(eventTimes);
            end
            
            % Initialize matrix for split events

            if(eventsExist)
                hasMergeEvents = ~cellfun('isempty',{tracks.mergeEvents})';
            end
            xMergeData = NaN(dLength, nTracks);
            yMergeData = NaN(dLength, nTracks);
            for i = find(hasMergeEvents)'
                eventTimes = false;
                eventTimes([tracks(i).mergeEvents tracks(i).mergeEvents+1]) = true;
                eventTimes = find(eventTimes)-1;
                dragtailWindow = [trackLengths(i) - displayLength(i) + 1 trackLengths(i)];
                eventTimes = eventTimes(eventTimes >= dragtailWindow(1) & eventTimes <= dragtailWindow(2));
                xMergeData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).xCoord(eventTimes);
                yMergeData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).yCoord(eventTimes);
            end
            
            % Plot tracks
            if isfield(tracks,'label') % If track is classified
                nColors = size(obj.Color,1);
                h = -ones(nColors,2);
                for iColor = 1:nColors
                    iTracks = mod([tracks.label]-1, nColors) +1 == iColor;
                    h(iColor,1)=plotFast(h(iColor,1),xData(:,iTracks),yData(:,iTracks),'Linestyle',obj.Linestyle,...
                        'Linewidth', obj.Linewidth, 'Color',obj.Color(iColor,:),varargin{:});
                    h(iColor,2)=plotFast(h(iColor,2),xGapData(:,iTracks),yGapData(:,iTracks),'Linestyle',obj.GapLinestyle,...
                        'Linewidth', obj.Linewidth, 'Color', obj.Color(iColor,:),varargin{:});
                end
            else
                % Plot links and gaps
                h=-ones(4,1);
                splitMarker = 'none';
                mergeMarker = 'none';
                if(obj.markMergeSplit)
                    splitMarker = obj.SplitMarker;
                    mergeMarker = obj.MergeMarker;
                end
                h(1) = plotFast(h(1),xData, yData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color',obj.Color,varargin{:});
                h(2) = plotFast(h(2),xGapData, yGapData, 'Linestyle', obj.GapLinestyle',...
                    'Linewidth', obj.Linewidth, 'Color',[1 1 1] - obj.Color, varargin{:});
                h(3) = plotFast(h(3),xSplitData, ySplitData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color', obj.SplitColor , 'Marker', splitMarker , varargin{:});
                h(4) = plotFast(h(4),xMergeData, yMergeData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color', obj.MergeColor, 'Marker', mergeMarker , varargin{:});
            end
            
            % Display track numbers if option is selected
            if obj.showLabel
                hlabels = -ones(nTracks,1);
                for i = find(~all(isnan(xData),1))
                    trackNr = num2str(tracks(i).number);
                    % Find last non-NaN coordinate
                    index = find(~isnan(xData(:,i)),1,'last');
                    if isfield(tracks,'label')
                        iColor = mod(tracks(i).label, nColors) + 1;
                        hlabels(i) = text(xData(index,i)+2, yData(index,i)+2, trackNr,...
                            'Color', obj.Color(iColor,:));
                    else
                        hlabels(i) = text(xData(index,i)+2, yData(index,i)+2, trackNr,...
                            'Color', obj.Color);
                    end
                end
                h = [h(:) ; hlabels ];
            end
            
            % Set tag
            set(h(ishandle(h)), 'Tag', tag, 'ButtonDownFcn', obj.ButtonDownFcn);
        end
        
        function updateDraw(obj, h, data)
            tag=get(h(1),'Tag');
            delete(h);
            obj.initDraw(data,tag);
            return;
            
        end
    end
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@(x)ischar(x) ||isvector(x);
            params(2).name='Linestyle';
            params(2).validator=@ischar;
            params(3).name='GapLinestyle';
            params(3).validator=@ischar;
            params(4).name='dragtailLength';
            params(4).validator=@isscalar;
            params(5).name='showLabel';
            params(5).validator=@isscalar;
            params(6).name='useDragtail';
            params(6).validator=@islogical;
            params(7).name='MergeColor';
            params(7).validator=@(x)ischar(x) ||isvector(x);
            params(8).name='SplitColor';
            params(8).validator=@(x)ischar(x) ||isvector(x);
            params(9).name='ButtonDownFcn';
            params(9).validator=@(x) isempty(x) || isa(x, 'function_handle');
            params(10).name='markMergeSplit';
            params(10).validator=@isscalar;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end
end
function h = plotFast(h,xData,yData,varargin)
%plotFast Uses the low-level plotting function line to plot matrices by
%separate lines into a single line with gaps while also trying to avoid
%creating a new object.
%
% INPUT
% h - handle to update if valid
% xData - same as given to plot, can be a matrix where columns are separate
% lines
% yData - same as given to plot, can be a matrix where columns are separate
% lines
% varargin - passes through extra parameters to the line function
%
% OUTPUT
%
% h - handle to the line object created or used
%

% Mark Kittisopikul, 2014/11/24

    xData(end+1,:) = NaN;
    yData(end+1,:) = NaN;
    if(ishandle(h))
        set(h,'XData',xData(:),'YData',yData(:),varargin{:});
    else
        h = line(xData(:),yData(:),varargin{:});
        if(isempty(h))
            h = -1;
        end
    end
end
