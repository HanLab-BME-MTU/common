classdef TracksDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Linestyle='-';
        GapLinestyle='--';
        Color='r';  
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
        function h=initDraw(obj,tracks,tag,varargin)
            nTracks = numel(tracks);
            h=-ones(nTracks,3);
                        
            % Get track length and display valid tracks
            trackLengths = cellfun(@numel,{tracks.tracksCoordAmpCG})/8;  
            validTracks = find(trackLengths>0);
            
            for i=validTracks
                xData= tracks(i).tracksCoordAmpCG(max(1,1+8*(trackLengths(i)-obj.dragtailLength)):8:end);
                yData= tracks(i).tracksCoordAmpCG(max(2,2+8*(trackLengths(i)-obj.dragtailLength)):8:end);
                % check gap is not exclusively composed of NaN's (e.g.
                % gaps with a small dragtail)
                if ~all(isnan(xData)) 
                    % Set color depending if track is classified or not
                    if isfield(tracks,'label')
                        color = obj.Color(mod(tracks(i).label,size(obj.Color,1))+1,:);
                    else
                        color =obj.Color;
                    end
                    
                    h(i,1)=plot(xData(~isnan(xData)),yData(~isnan(yData)),...
                        'Linestyle',obj.GapLinestyle','Color',color,...
                        varargin{:});
                    h(i,2)=plot(xData,yData,'Linestyle',obj.Linestyle,...
                        'Color',color,varargin{:});
                    if obj.showLabel
                        h(i,3) = text(xData(find(~isnan(xData),1,'last'))+2,...
                            yData(find(~isnan(yData),1,'last'))+2,num2str(i),...
                            'Color',color);
                    end
                end
            end
            set(h(ishandle(h)),'Tag',tag); 
           
        end

        function updateDraw(obj,allh,data)
            tag=get(allh(1),'Tag');
            delete(allh);
            obj.initDraw(data,tag);
            return;
            nTracks = numel(data);

            h=findobj(allh,'Type','line');
            delete(h(2*nTracks+1:end));
            h(2*nTracks+1:end)=[];
            hlinks=findobj(h,'LineStyle',obj.Linestyle);
            hgaps=findobj(h,'LineStyle',obj.GapLinestyle);
            
            % Update existing windows
            for i=1:min(numel(hlinks),nTracks) 
                xData= data.x{i}(max(1,end-obj.dragtailLength):end);
                yData= data.y{i}(max(1,end-obj.dragtailLength):end);
                set(hgaps(i),'Xdata',xData(~isnan(xData)),'YData',yData(~isnan(yData)));
                set(hlinks(i),'Xdata',xData,'YData',yData);
            end
            for i=min(numel(hlinks),nTracks)+1:nTracks
                xData= data.x{i}(max(1,end-obj.dragtailLength):end);
                yData= data.y{i}(max(1,end-obj.dragtailLength):end);
                hgaps(i)=plot(xData(~isnan(xData)),yData(~isnan(yData)),...
                    'Linestyle',obj.GapLinestyle','Color',obj.Color);
                hlinks(i)=plot(xData,yData,...
                    'Linestyle',obj.Linestyle,'Color',obj.Color);
            end
            set([hlinks hgaps],'Tag',tag);

            
            if isfield(data,'label') && obj.showLabel
                ht=findobj(allh,'Type','text');
                delete(ht(nTracks+1:end));
                ht(nTracks+1:end)=[];
                for i=1:min(numel(ht),nTracks)
                    set(ht(i),'Position',[data.x{i}(end),data.y{i}(end)],...
                        'String',data.label(i));
                end
                for i=min(numel(ht),nTracks)+1:nTracks
                    ht(i) = text(data.x{i}(end),data.y{i}(end),num2str(data.label(i)));
                end
                set(ht,'Tag',tag); 
            end
           
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
        end

        function f=getDataValidator() 
            f=@isstruct;
        end
    end    
end