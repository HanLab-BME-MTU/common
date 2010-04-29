classdef SegmentationProcess < Process
% A concrete process for mask process info
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
       maskPaths_ 
       funName_
       funParams_
    end
    
    methods (Access = public)
        function obj = SegmentationProcess (owner, maskPaths, ...
                funName, funParams)
           % Construntor of class MaskProcess
           if nargin == 0
              super_args = {};
           else
               super_args{1} = owner;
               super_args{2} = 'SegmentationProcess'; % ??? process name/ID ???
           end
           % Call the supercalss constructor with empty cell array (no
           % argument) if nargin == 0
           obj = obj@Process(super_args{:});
           if nargin > 0
              obj.maskPaths_ = maskPaths;
              obj.funName_ = funName;
              obj.funParams_ = funParams;
           end
        end
        function sanityCheck(obj)
            % Sanity Check

            if length(obj.maskPath_) ~= length(obj.owner_.channelPath_)
                
                disp('User-defined warning: The number of mask channels is different from input channel');
            end
            % Check mask path for each channel
            % ... ...
            
        end
    end
end