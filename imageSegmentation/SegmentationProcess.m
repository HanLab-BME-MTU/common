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
               super_args{2} = 'Segmentation'; 
           end
           % Call the supercalss constructor with empty cell array (no
           % argument) if nargin == 0
           obj = obj@Process(super_args{:});
           if nargin > 0
              obj.maskPaths_ = maskPaths;
              obj.funName_ = funName;
              if isnumeric (funParams)
                obj.funParams_ = funParams;
              else
                obj.funParams_ = str2double(funParams);
              end
           end
        end
        function setPara(obj, para)
            % Reset process' parameters
            obj.funParams_ = para;
        end
        function sanityCheck(obj) % throw exception
            % Sanity Check
            if obj.funParams_ < 0
                error('lccb:set:fatal','Fatal problem is detected in pamameter set up.\n\n');
            end
            % Check mask path for each channel
            % ... ...
        end
    end
    methods (Static)
        function text = getHelp(obj)
           text = 'This is specific help text of segmentation process'; 
        end
    end
end