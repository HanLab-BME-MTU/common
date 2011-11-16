classdef DetectionProcess < ImageAnalysisProcess
    % A class definition for a generic detection process.
    %
    % Chuangang Ren, 11/2010
    
    properties(SetAccess = protected, GetAccess = public)
        
        channelIndex_ % The index of channel to process
        overwrite_ = 0; % If overwrite the original result MAT file
    end
    
    methods(Access = public)
        
        function obj = DetectionProcess(owner, name, funName, channelIndex, funParams )
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;
            end
            
            if nargin > 3
                obj.channelIndex_ = channelIndex;
            end
            
            if nargin > 4
                obj.funParams_ = funParams;
            end
            
            
        end
        
        
        % Set overwrite
        function setOverwrite (obj, i)
            obj.overwrite_ = i;
        end
        
        
        function setChannelIndex(obj, index)
            % Set channel index
            if any(index > length(obj.owner_.channels_))
                error ('User-defined: channel index is larger than the number of channels.')
            end
            if ~isequal(obj.channelIndex_,index)
                obj.channelIndex_ = index;
                obj.procChanged_=true;
            end
        end
        
    
        
    end
    methods(Static)
        function name = getName()
            name = 'Detection';
        end
        function h = GUI()
            h= @detectionProcessGUI;
        end
    end

end