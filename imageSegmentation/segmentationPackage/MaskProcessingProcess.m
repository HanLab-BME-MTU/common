classdef MaskProcessingProcess < SegmentationProcess
    %Generic class definition for processes which use / post-process / edit
    %masks which have already been created, but do not themselves directly
    %create masks.
    %
    
    properties (SetAccess = private,GetAccess = public)
    
       inMaskPaths_ %Mask directories which were used as input.              
        
    end
    
    methods(Access = public)
        
        function obj = MaskProcessingProcess(owner,name,funName, funParams,...
                                             inMaskPaths,outMaskPaths)
           % Constructor of class MaskProcessingProcess
           if nargin == 0
              super_args = {};
           else
               super_args{1} = owner;
               super_args{2} = name;                
           end
           if nargin > 2
               super_args{3} = funName;
           end
           if nargin > 3
               super_args{4} = funParams;
           end
           if nargin > 5
               super_args{5} = outMaskPaths;
           end
           
           % Call the superclass constructor - these values are private
           obj = obj@SegmentationProcess(super_args{:});
           
           if nargin > 4               
              if ~isempty(inMaskPaths) && numel(inMaskPaths) ...
                      ~= numel(owner.channels_) || ~iscell(inMaskPaths)
                 error('lccb:set:fatal','Mask paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end
              obj.inMaskPaths_ = inMaskPaths;              
           else
               obj.inMaskPaths_ = cell(1,numel(owner.channels_));               
           end
        end
%         function sanityCheck(obj) % throws exception
%             % Sanity Check
%             % 1. Check that non-empty mask paths exist
%             % 2. Check number of mask = number of raw images
%             
%             for i = 1:length(obj.inMaskPaths_)
%                 if ~isempty(obj.inMaskPaths_{i})
%                     if ~exist(obj.inMaskPaths_{i}, 'dir')
%                         error('lccb:set:fatal','Cannot find mask paths %s.\n\n', ...
%                                 obj.inMaskPaths_{i});
%                     end                                        
%                     
%                     % Number of mask image files
%                     maskFileNames = imDir(obj.inMaskPaths_{i}, true);
%                     
%                     if length(maskFileNames) ~= obj.owner_.nFrames_;
%                         error('lccb:set:fatal', 'Invalid mask input directory: The number of masks in %s is inconsistent with the number of input images in %s.\n\n',...
%                             obj.inMaskPaths_{i}, obj.owner_.getChannelPaths(i);
%                     end
%                 end
%             end
%         end        
        function setInMaskPath(obj,iChan,maskPaths)           
            if all(obj.checkChanNum(iChan))
                nChan = numel(iChan);
                if ~iscell(maskPaths)
                    maskPaths = {maskPaths};
                end
                for j = 1:nChan
                    if exist(maskPaths{j},'dir') && ....
                            length(imDir(maskPaths{j})) == obj.owner_.nFrames_;
                        obj.inMaskPaths_{iChan(j)} = maskPaths{j};
                    else
                        error(['lccb:set:fatal','Invalid input mask path for channel ' num2str(iChan(j))]);
                    end
                end
            else
                error('lccb:set:fatal','Invalid mask channel number for mask path!\n\n'); 
            end
        end
        function fileNames = getInMaskFileNames(obj,iChan)
            if all(obj.checkChanNum(iChan))
                fileNames = cellfun(@(x)(imDir(x)),obj.inMaskPaths_(iChan),'UniformOutput',false);
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

end