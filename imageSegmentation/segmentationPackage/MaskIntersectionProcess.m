classdef MaskIntersectionProcess < MaskProcessingProcess
    % A concrete process to create mask interesections
    
    
    methods(Access = public)        
        function obj = MaskIntersectionProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = MaskIntersectionProcess.getName;
                super_args{3} = @intersectMovieMasks;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%                                            
                    funParams.ChannelIndex = 1:nChan; %Default is to transform masks for all channels
                    funParams.SegProcessIndex = []; %No default...
                    funParams.OutputDirectory =  [outputDir filesep 'intersected_masks'];      
                    funParams.BatchMode = false;                                              
                    
                end
                
                super_args{4} = funParams;                    
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
            
        end  

        
        %Checks if a particular channel has masks
        function maskStatus = checkChannelOutput(obj)
            maskStatus = checkChannelOutput@MaskProcessingProcess(obj,1);
        end
        
        function mask = loadChannelOutput(obj,iFrame,varargin)
            mask = loadChannelOutput@MaskProcessingProcess(obj,1,iFrame,varargin{:});
        end
        
        function h=draw(obj,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iFrame',@isnumeric);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,iFrame,varargin{:})
			
            data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput}));
            catch ME
                obj.displayMethod_{iOutput}=...
                    outputList(iOutput).defaultDisplayMethod();
            end
            
            % Delegate to the corresponding method
            tag = [obj.getName '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end
        
        function output = getDrawableOutput(obj)
            output=getDrawableOutput@MaskProcessingProcess(obj);
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@LineDisplay;
        end
        
        
        
        
    end
    methods(Static)
        function name =getName()
            name = 'Mask Intersection';
        end
    end
end
    