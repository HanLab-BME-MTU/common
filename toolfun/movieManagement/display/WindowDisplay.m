classdef WindowDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='r';  
        FaceAlpha=.2;
        showNum=0;
    end
    methods
        function obj=VectorFieldDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            for j = 1:numel(data)
                for k = 1:numel(data{j})
                    if ~isempty(data{j}{k})
                        currWin = [data{j}{k}{:}];
                        if ~isempty(currWin)
                            
                            h(j,k)=patch(currWin(1,:),currWin(2,:),obj.Color,...
                                'FaceAlpha',obj.FaceAlpha);
                            
%                             if showNum && mod(j,showNum)==0 && mod(k,showNum) == 0
%                                 text(currWin(1,1),currWin(2,1),[num2str(j) ',' num2str(k)])
%                             end
                        end
                    end
                end
            end
%             h=quiver(data(:, 2),data(:, 1),obj.scale*(data(:, 4)-data(:, 2)),...
%                 obj.scale*(data(:, 3)-data(:, 1)),0,'-','Autoscale',autoscale,...
%                 'Color',obj.Color,varargin{:});
            set(nonzeros(h(:)),'Tag',tag);
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.FaceAlpha=ip.Results.FaceAlpha;
        end
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            delete(h);
            h=[];
            for j = 1:numel(data)
                for k = 1:numel(data{j})
                    if ~isempty(data{j}{k})
                        currWin = [data{j}{k}{:}];
                        if ~isempty(currWin)
                            
                            h(j,k)=patch(currWin(1,:),currWin(2,:),obj.Color,...
                                'FaceAlpha',obj.FaceAlpha);
                            
%                             if showNum && mod(j,showNum)==0 && mod(k,showNum) == 0
%                                 text(currWin(1,1),currWin(2,1),[num2str(j) ',' num2str(k)])
%                             end
                        end
                    end
                end
            end
            set(nonzeros(h(:)),'Tag',tag);
%             set(h,'XData',data(:,2),'YData',data(:,1),...
%                 'UData',obj.scale*(data(:, 4)-data(:, 2)),...
%                 'VData',obj.scale*(data(:, 3)-data(:, 1)))
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
            ip.addParamValue('FaceAlpha',obj.FaceAlpha,@isscalar);
            ip.addParamValue('showNum',obj.showNum,@isscalar);
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@iscell;
        end
    end    
end