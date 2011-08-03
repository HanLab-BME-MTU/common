classdef CircleDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Color='r';        
    end
    methods
        function obj=LineDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            % Calculate sin and cos vectors to display circles
            ang=0:pi/4:2*pi;
            sn_small=sin(ang);
            cs_small=cos(ang);
            ang=0:pi/10:2*pi;
            sn_med=sin(ang);
            cs_med=cos(ang);
            ang=0:pi/100:2*pi;
            sn_large=sin(ang);
            cs_large=cos(ang);
            
            smallInd = data(:,3)<=2;
            largeInd = data(:,3)>10;
            medInd=~smallInd&&~largeInd;
            h(1)=plot(data(smallInd,2)+data(smallInd,3)*sn_small,...
                data(smallInd,1)+data(smallInd,3)*cs_small,'.');
            h(2)=plot(data(medInd,2)+data(medInd,3)*sn_med,...
                data(medInd,1)+data(medInd,3)*cs_med,'r.');
            h(3)=plot(data(largeInd,2)+data(largeInd,3)*sn_large,...
                data(largeInd,1)+data(largeInd,3)*cs_large,'.');
            set(h,'Tag',tag,'Color',obj.Color);
        end
        function updateDraw(obj,h,data)
            ang=0:pi/4:2*pi;
            sn_small=sin(ang);
            cs_small=cos(ang);
            ang=0:pi/10:2*pi;
            sn_med=sin(ang);
            cs_med=cos(ang);
            ang=0:pi/100:2*pi;
            sn_large=sin(ang);
            cs_large=cos(ang);
            
            smallInd = data(:,3)<=2;
            largeInd = data(:,3)>10;
            medInd=~smallInd&&~largeInd;        
            set(h(1),'XData',data(smallInd,2)+data(smallInd,3)*sn_small,...
                'YData',data(smallInd,1)+data(smallInd,3)*cs_small);
            set(h(2),'XData',data(medInd,2)+data(medInd,3)*sn_med,...
                'YData',data(medInd,1)+data(medInd,3)*cs_med);
            set(h(3),'XData',data(largeInd,2)+data(largeInd,3)*sn_large,...
                'YData',data(largeInd,1)+data(largeInd,3)*cs_large);
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@ischar);
        end 
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
        end
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end