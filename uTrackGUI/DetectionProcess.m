classdef DetectionProcess < Process
% A class definition for a generic detection process.
%
% Chuangang Ren, 11/2010
 
properties(SetAccess = protected, GetAccess = public)
    
    outParams_ % All output data or path
    channelIndex_ % The index of channel to process
end

methods(Access = public)
    
    function obj = DetectionProcess(owner, name, funName, channelIndex, funParams )
       
        if nargin == 0
            super_args = {};
        else
            super_args{1} = owner;
            super_args{2} = name;
        end
        
        obj = obj@Process(super_args{:});
        
        if nargin > 2
           obj.funName_ = funName; 
        end
        
        if nargin > 3
           obj.channelIndex_ = channelIndex; 
        end
        
        if nargin > 4
            obj.funParams_ = funParams;
        end
        
        obj.outParams_ = cell(1, length(owner.channels_));

    end
    
    function sanityCheck(obj) 
    end
    
    function setOutPara(obj, para)
        % Reset process' parameters
        obj.outParams_ = para;
    end
    
    function setChannelIndex(obj, index)
        % Set channel index
        if any(index > length(obj.owner_.channels_))
           error ('User-defined: channel index is larger than the number of channels.') 
        end
        obj.channelIndex_ = index;
    end
    
end


end