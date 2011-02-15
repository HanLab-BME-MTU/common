classdef ImageAnalysisProcess < Process
%A class definition for a generic image analysis process. That is, a
%process which takes in images and produces some other (non-image) form 
%of output. This generic class expects there to be one file per frame, per
%channel with each channel in a separate directory.
%
%
% Hunter Elliott, 7/2010
%
    properties (SetAccess = protected, GetAccess = public)
        
        outFilePaths_ %Location of result file(s) for each channel.
        inImagePaths_ %Location of images which were processed. May or may not be the original, raw images.        
        
    end
    
    methods (Access = public)
        
        function obj = ImageAnalysisProcess(owner,name,funName,funParams,...
                                            inImagePaths,outFilePaths)
                                          
            if nargin == 0;
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
               obj.funParams_ = funParams;              
            end
            
            if nargin > 4
              if ~isempty(inImagePaths) && numel(inImagePaths) ...
                      ~= numel(owner.channels_) || ~iscell(inImagePaths)
                 error('lccb:set:fatal','Input image paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end
              obj.inImagePaths_ = inImagePaths;              
            else
                %Default is to use raw images as input.
                obj.inImagePaths_ = owner.getChannelPaths;               
            end                        
            if nargin > 5               
              if ~isempty(outFilePaths) && numel(outFilePaths) ...
                      ~= numel(owner.channels_) || ~iscell(outFilePaths)
                 error('lccb:set:fatal','Output File paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end
              obj.outFilePaths_ = outFilePaths;              
            else
                obj.outFilePaths_ = cell(1,numel(owner.channels_));               
            end
            
        end
        
        function setOutFilePath(obj,chanNum,filePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for file path!\n\n'); 
            end
            
            if ~iscell(filePath)
                filePath = {filePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(filePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
               if ~exist(filePath{j},'dir') && ~exist(filePath{j},'file')
                   error('lccb:set:fatal',...
                       ['The directory specified for output for channel ' ...
                       num2str(chanNum(j)) ' is invalid!']) 
               else
                   obj.outFilePaths_{chanNum(j)} = filePath{j};                
               end
            end
            
            
        end
        function setInImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n'); 
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
               if ~exist(imagePath{j},'dir')
                   error('lccb:set:fatal',...
                       ['The directory specified for channel ' ...
                       num2str(chanNum(j)) ' is invalid!']) 
               
               else
                   if isempty(imDir(imagePath{j}))
                       error('lccb:set:fatal',...
                       ['The directory specified for channel ' ...
                       num2str(chanNum(j)) ' does not contain any images!!']) 
                   else                       
                       obj.inImagePaths_{chanNum(j)} = imagePath{j};                
                   end
               end
            end                        
        end
        
        function fileNames = getInImageFileNames(obj,iChan)
            if obj.checkChanNum(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.inImagePaths_(iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)                    
                    error('Incorrect number of images found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end    
            
            
        end
        
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channels have valid output files
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan)
               iChan = 1:nChanTot;
           end
           %Makes sure there's at least one .mat file in the speified
           %directory
           OK =  arrayfun(@(x)(x <= nChanTot && ...
                             x > 0 && isequal(round(x),x) && ...
                             exist(obj.outFilePaths_{x},'dir') && ...
                             ~isempty(dir([obj.outFilePaths_{x} filesep '*.mat']))),iChan);
        end
        
        function figHan = showResult(obj)
            %There is no generic showResult method for this class
            figHan = [];            
            
        end
               
            
        function sanityCheck(obj)
            
        end
        
        
        
        
    end
    
end