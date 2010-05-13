classdef SegmentationProcess < Process
% A concrete process for mask process info
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
       maskPaths_

    end
    
    methods (Access = public)
        function obj = SegmentationProcess (owner,funName, funParams,...
                        maskPaths,outputDirectory)
           % Constructor of class MaskProcess
           if nargin == 0
              super_args = {};
           else
               super_args{1} = owner;
               super_args{2} = 'Segmentation'; 
               super_args{3} = outputDirectory;
           end
           % Call the superclass constructor with empty cell array (no
           % argument) if nargin == 0
           obj = obj@Process(super_args{:});
           if nargin > 0
              if ~isempty(maskPaths) && numel(maskPaths) ...
                      ~= numel(owner.channelPath_) || ~iscell(maskPaths)
                 error('lccb:set:fatal','Mask paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end
              obj.maskPaths_ = maskPaths;
              obj.funName_ = funName;
              obj.funParams_ = funParams;      
              
           end
        end
        function sanityCheck(obj) % throw exception
            % Sanity Check
            if obj.funParams_ < 0
                error('lccb:set:fatal','Function parameter should be larger than 0.\n\n');
            end
            % Check mask path for each channel
            % ... ...
        end
        function setMaskPath(obj,chanNum,maskPath)           
            if isnumeric(chanNum) && chanNum > 0 && chanNum < numel(owner.channelPaths_)
                obj.maskPaths_{chanNum} = maskPath;
            else
                error('lccb:set:fatal','Invalid mask channel number for mask path!\n\n'); 
            end
        end
    end
    methods (Static)
        function text = getHelp(obj)
           text = 'This process will create masks for the selected movie channels. These masks will be saved to a directory specified by the user as binary .tif files.'; 
        end
    end
end