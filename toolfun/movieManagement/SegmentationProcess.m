classdef SegmentationProcess < Process
% A concrete process for mask process info
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
       maskPaths_

    end
    
    methods (Access = public)
        function obj = SegmentationProcess (owner,name,funName, funParams,...
                        maskPaths)
           % Constructor of class SegmentationProcess
           if nargin == 0
              super_args = {};
           else
               super_args{1} = owner;
               super_args{2} = name;                
           end
           % Call the superclass constructor - these values are private
           obj = obj@Process(super_args{:});

           if nargin > 2
               obj.funName_ = funName;                              
           end
           if nargin > 3
              obj.funParams_ = funParams;              
           end
           if nargin > 4               
              if ~isempty(maskPaths) && numel(maskPaths) ...
                      ~= numel(owner.channelPath_) || ~iscell(maskPaths)
                 error('lccb:set:fatal','Mask paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end
              obj.maskPaths_ = maskPaths;              
           else
               obj.maskPaths_ = cell(1,numel(owner.channelPath_));               
           end
        end
        function sanityCheck(obj) % throws exception
            % Sanity Check
            % 1. Check that non-empty mask paths exist
            % 2. Check number of mask = number of raw images
            
            for i = 1: length(obj.maskPaths_)
                if ~isempty(obj.maskPaths_{i})
                    if ~exist(obj.maskPaths_{i}, 'dir')
                        error('lccb:set:fatal','Cannot find mask paths %s.\n\n', ...
                                obj.maskPaths_{i});
                    end
                    
                    % Number of raw image files
                    fileNames = imDir(obj.owner_.channelPath_{i}, true);
                    
                    % Number of mask image files
                    maskFileNames = imDir(obj.maskPaths_{i}, true);
                    
                    if length(fileNames) ~= length(maskFileNames)
                        error('lccb:set:fatal', 'The number of masks in %s is inconsistent with the number of input images in %s.\n\n',...
                            obj.maskPaths_{i}, obj.owner_.channelPath_{i});
                    end
                end
            end
        end
        function setMaskPath(obj,chanNum,maskPath)           
            if isnumeric(chanNum) && chanNum > 0 && ...
                    chanNum <= numel(obj.owner_.channelPath_)
                obj.maskPaths_{chanNum} = maskPath;
            else
                error('lccb:set:fatal','Invalid mask channel number for mask path!\n\n'); 
            end
        end
        function fileNames = getMaskFileNames(obj,iChan)
            if isnumeric(iChan) && min(iChan)>0 && max(iChan) <= ...
                    numel(obj.owner_.channelPath_) && isequal(round(iChan),iChan)                
                fileNames = cellfun(@(x)(imDir(x)),obj.maskPaths_(iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)                    
                    error('Incorrect number of masks found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end    
            
            
        end
    end
    methods (Static)
        function text = getHelp(obj)
           text = 'This process will create masks for the selected movie channels. These masks will be saved to a directory specified by the user as binary .tif files.'; 
        end
    end
end