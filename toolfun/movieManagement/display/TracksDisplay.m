classdef TracksDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Linestyle='-';
        Linewidth=1;
        GapLinestyle='--';
        Color='r';
        MergeColor = 'y';
        SplitColor = 'g';
        useDragtail=true;
        dragtailLength=10;
        showLabel=false;
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
            trackLengths = cellfun(@numel,{tracks.xCoord});
            validTracks = find(trackLengths>0);
            tracks = tracks(validTracks);
            trackLengths = trackLengths(validTracks);
            
            nTracks = numel(validTracks);
            h=-ones(nTracks,3);
            
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
            for i = 1 : nTracks
                displayLength = trackLengths(i) - max(0,trackLengths(i)-dLength);
                xData(1:displayLength, i) = tracks(i).xCoord(end-displayLength+1:end);
                yData(1:displayLength, i) = tracks(i).yCoord(end-displayLength+1:end);
            end
            
            % Initialize matrix for gaps
            xGapData = NaN(size(xData));
            yGapData = NaN(size(xData));
            
            % Label gaps: series of NaNs not-connected to the border
            I = isnan(xData);
            I = [I; zeros(size(I))];
            I = reshape(I, size(I,1)/2, size(I,2)*2);
            I = [zeros(size(I,1), 1) I];
            I = bwlabel(imclearborder(I));
            I = I(:, 2:2:end);
            
            % Fill gaps x and y data
            for i = unique(nonzeros(I))'
                iFirst = find(I == i, 1, 'first')-1;
                iLast = find(I == i, 1, 'last')+1;
                xGapData(iFirst:iLast) = linspace(xData(iFirst), xData(iLast), iLast - iFirst +1);
                yGapData(iFirst:iLast) = linspace(yData(iFirst), yData(iLast), iLast - iFirst +1);
            end
            
            % Initialize matrix for split events
            hasSplitEvents = arrayfun(@(x) ~isempty(find(x.events=='s',1)), tracks);
            xSplitData = NaN(dLength, nTracks);
            ySplitData = NaN(dLength, nTracks);
            for i = find(hasSplitEvents)'
                eventTimes = find(tracks(i).events == 's');
                eventTimes = unique([eventTimes, eventTimes+1]);
                dragtailWindow = trackLengths(i)-displayLength+1:trackLengths(i);
                eventTimes = intersect(eventTimes, dragtailWindow);
                xSplitData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).xCoord(eventTimes);
                ySplitData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).yCoord(eventTimes);
            end
            
            % Initialize matrix for split events
            hasMergeEvents = arrayfun(@(x) ~isempty(find(x.events=='m',1)), tracks);
            xMergeData = NaN(dLength, nTracks);
            yMergeData = NaN(dLength, nTracks);
            for i = find(hasMergeEvents)'
                eventTimes = find(tracks(i).events == 'm');
                eventTimes = unique([eventTimes, eventTimes - 1]);
                dragtailWindow = trackLengths(i)-displayLength+1:trackLengths(i);
                eventTimes = intersect(eventTimes, dragtailWindow);
                xMergeData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).xCoord(eventTimes);
                yMergeData(eventTimes - dragtailWindow(1) +1, i) = tracks(i).yCoord(eventTimes);
            end
            
            % Plot tracks
            if isfield(tracks,'label') % If track is classified
                nColors = size(obj.Color,1);
                for iColor = 1:nColors
                    iTracks = mod([tracks.label]-1, nColors) +1 == iColor;
                    h(iTracks,1)=plot(xData(:,iTracks),yData(:,iTracks),'Linestyle',obj.Linestyle,...
                        'Linewidth', obj.Linewidth, 'Color',obj.Color(iColor,:),varargin{:});
                    h(iTracks,2)=plot(xGapData(:,iTracks),yGapData(:,iTracks),'Linestyle',obj.GapLinestyle,...
                        'Linewidth', obj.Linewidth, 'Color', obj.Color(iColor,:),varargin{:});
                end
            else
                % Plot links and gaps
                h(:,1) = plot(xData, yData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color',obj.Color,varargin{:});
                h(:,2) = plot(xGapData, yGapData, 'Linestyle', obj.GapLinestyle',...
                    'Linewidth', obj.Linewidth, 'Color',[1 1 1] - obj.Color, varargin{:});
                h(:,3) = plot(xSplitData, ySplitData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color','y', varargin{:});
                h(:,4) = plot(xMergeData, yMergeData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color','g', varargin{:});
            end
            
            % Display track numbers if option is selected
            if obj.showLabel
                for i = find(~all(isnan(xData),1))
                    trackNr = num2str(tracks(i).number);
                    % Find last non-NaN coordinate
                    index = find(~isnan(xData(:,i)),1,'last');
                    if isfield(tracks,'label')
                        iColor = mod(tracks(i).label, nColors) + 1;
                        h(i,3) = text(xData(index,i)+2, yData(index,i)+2, trackNr,...
                            'Color', obj.Color(iColor,:));
                    else
                        h(i,3) = text(xData(index,i)+2, yData(index,i)+2, trackNr,...
                            'Color', obj.Color);
                    end
                end
            end
            
            % Set tag
            set(h(ishandle(h)), 'Tag', tag);
            
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
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end
end
