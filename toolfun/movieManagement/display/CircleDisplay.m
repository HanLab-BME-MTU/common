classdef CircleDisplay < MovieDataDisplay
    %Concrete display class for displaying errors as circles
    % Adapated from fsmVectorAnalysis
    properties
        Color='r';        
    end
    methods
        function obj=CircleDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            % Convert data into circles
            [x,y] = genCircles(data);       
       
            % Plot circles and set tag
            h=plot(x', y','.','Color',obj.Color);
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            % Retrieve the tag
            tag=get(h(1),'Tag');
            
            % Convert data into circles
            [x,y] = genCircles(data);            
            
            % Delete extra plots handles
            delete(h(size(x,1)+1:end));
            h(size(x,1)+1:end)=[];
            
            % Update existing handles
            for i=1:min(numel(h),size(x,1))
                set(h(i),'Xdata',x(i,:),'YData',y(i,:));
            end
            
            % Plot additional circles
            addIndx= min(numel(h),size(x,1))+1:size(x,1);
            h(addIndx)=plot(x(addIndx,:)', y(addIndx,:)','.','Color',obj.Color);
            
            % Set tag
            set(h,'Tag',tag);
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
            f=@(x)isnumeric(x) & size(x,2)>=3 ;
        end
    end    
end

function [x,y] = genCircles(data)
ang=0:pi/100:2*pi;
sn=sin(ang);
cs=cos(ang);

x=repmat(data(:,2),1,size(sn,2))+data(:,3)*sn;
y=repmat(data(:,1),1,size(cs,2))+data(:,3)*cs;

end