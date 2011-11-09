classdef TracksDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Linestyle='-';
        Color='r';  
        dragtailLength=10;
        showLabel=true;
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
        function h=initDraw(obj,data,tag,varargin)
            nTracks = numel(data.x);
            h=-ones(nTracks,1);
            for i=1:nTracks
                h(i)=plot(data.x{i}(max(1,end-obj.dragtailLength):end),...
                    data.y{i}(max(1,end-obj.dragtailLength):end),varargin{:});
            end
            set(h,'Tag',tag,'LineStyle',obj.Linestyle,'Color',obj.Color); 
            
            if isfield(data,'label')  && obj.showLabel
                ht=-ones(nTracks,1);
                for i=1:nTracks
                    ht(i) = text(data.x{i}(end),data.y{i}(end),num2str(data.label(i)));
                end
                set(ht,'Tag',tag);

            end
        end

        function updateDraw(obj,allh,data)
            tag=get(allh(1),'Tag');
            nTracks = numel(data.x);

            h=allh(strcmp(get(allh,'Type'),'line'));
            ht=allh(strcmp(get(allh,'Type'),'text'));
            delete(h(nTracks+1:end));
            h(nTracks+1:end)=[];


            % Update existing windows
            for i=1:min(numel(h),nTracks)
                set(h(i),'Xdata',data.x{i}(max(1,end-obj.dragtailLength):end),...
                    'YData',data.y{i}(max(1,end-obj.dragtailLength):end));
            end
            for i=min(numel(h),nTracks)+1:nTracks
                h(i)=plot(data.x{i}(max(1,end-obj.dragtailLength):end),...
                    data.y{i}(max(1,end-obj.dragtailLength):end));
            end
            set(h,'Tag',tag,'Linestyle',obj.Linestyle,'Color',obj.Color);

            if isfield(data,'label') && obj.showLabel
                delete(ht(nTracks+1:end));
                ht(nTracks+1:end)=[];
                for i=1:min(numel(ht),nTracks)
                    set(ht(i),'Position',[data.x{i}(end),data.y{i}(end)],...
                        'String',data.label(i));
                end
                for i=min(numel(ht),nTracks)+1:nTracks
                    ht(i) = text(data.x{i}(end),data.y{i}(end),num2str(data.label(i)));
                end
            end
            set(ht,'Tag',tag); 
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.Linestyle=ip.Results.Linestyle;
            obj.dragtailLength=ip.Results.dragtailLength;
            obj.showLabel=ip.Results.showLabel;
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
            ip.addParamValue('Linestyle',obj.Linestyle,@ischar);
            ip.addParamValue('dragtailLength',obj.dragtailLength,@isscalar);
            ip.addParamValue('showLabel',obj.showLabel,@isscalar);
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isstruct;
        end
    end    
end