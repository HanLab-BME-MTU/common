classdef WindowsDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='r';  
        FaceAlpha=.2;
        showNum=5;
    end
    methods
        function obj=WindowsDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            
            [W,j,k] = prepareData(data);
            
            % Plot the windows as patches
            hw=cellfun(@(w)patch(w(1,:),w(2,:),obj.Color,...
                'FaceAlpha',obj.FaceAlpha),W);
            
            % Find indices of thes windows numbers to plot
            textInd = find(mod(j,obj.showNum)==0 & mod(k,obj.showNum)==0);            
            ht = arrayfun(@(i)text(W{i}(1,1),W{i}(2,1),...
                [num2str(j(i)) ',' num2str(k(i))]),textInd);

            h=nonzeros([hw ht]);
            set(h,'Tag',tag); 
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.FaceAlpha=ip.Results.FaceAlpha;
        end
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            
            [W,j,k] = prepareData(data);
            
            hw=h(strcmp(get(h,'Type'),'patch'));
            ht=h(strcmp(get(h,'Type'),'text'));
            
            % Delete extra windows
            delete(hw(numel(W)+1:end));
            hw(numel(W)+1:end)=[];
            
            % Update existing windows
            for i=1:min(numel(hw),numel(W))
                set(hw(i),'Xdata',W{i}(1,:),'YData',W{i}(2,:));
            end
            
            % Plot additional windows
            addWindows= min(numel(hw),numel(W))+1:numel(W);
            hw(addWindows)=cellfun(@(w)patch(w(1,:),w(2,:),obj.Color,...
                'FaceAlpha',obj.FaceAlpha),W(addWindows));
            
            %Retrieve text index
            textInd = find(mod(j,obj.showNum)==0 & mod(k,obj.showNum)==0);
            % Delete extra text objects
            delete(ht(numel(textInd)+1:end));
            ht(numel(textInd)+1:end)=[];
             
            % Update existing text objects
            for i=1:min(numel(ht),numel(textInd))
                wInd=textInd(i);
                set(ht(i),'Position',[W{wInd}(1,1) W{wInd}(2,1) 0],...
                    'String',[num2str(j(wInd)) ',' num2str(k(wInd))]);
            end
            
            % Plot additional text objects
            addText= min(numel(ht),numel(textInd))+1:numel(textInd);
            ht(addText) = arrayfun(@(i)text(W{i}(1,1),W{i}(2,1),...
                [num2str(j(i)) ',' num2str(k(i))]),textInd(addText));
            
            h=nonzeros([hw(:); ht(:)]);
            set(h,'Tag',tag); 

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

function [W,j,k] = prepareData(data)
% Concatenate the windows
W=cellfun(@(x) cellfun(@(y) [y{:}],x,'Unif',false),...
    data,'Unif',false);
W=horzcat(W{:});

% Create corresponding indices for labeling
j=arrayfun(@(x) repmat(x,1,numel(data{x})),1:numel(data),'Unif',false);
j=horzcat(j{:});
k=cellfun(@(x) 1:numel(x),data,'Unif',false);
k=horzcat(k{:});
end