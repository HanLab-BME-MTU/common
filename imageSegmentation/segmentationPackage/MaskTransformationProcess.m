classdef MaskTransformationProcess < MaskProcessingProcess
    %Class definition for post processing using refineMovieMasks.m
    
    
    methods(Access = public)        
        function obj = MaskTransformationProcess(owner,outputDir,funParams)
            
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = MaskTransformationProcess.getName;
                super_args{3} = @transformMovieMasks;                               
                
                if nargin < 3 || isempty(funParams)                                       
                    
                    %----Defaults----%                                            
                    funParams.ChannelIndex = 1:nChan; %Default is to transform masks for all channels
                    funParams.SegProcessIndex = []; %No default...
                    funParams.OutputDirectory = [outputDir filesep 'transformed_masks'];      
                    funParams.TransformFilePaths = cell(1,numel(owner.channels_));%No default...
                    funParams.BatchMode = false;                                              
                    
                end
                
                super_args{4} = funParams;                    
            end
            
            obj = obj@MaskProcessingProcess(super_args{:});
            
        end  
        function setTransformFilePath(obj,iChan,transformPath)                                   
        
            %Make sure the specified channels are valid
            if ~obj.checkChanNum(iChan)
                error('The channel indices specified for transform files are invalid!')
            end
                
            %If only one transform path was input, convert to a cell
            if ~iscell(transformPath)
                transformPath = {transformPath};
            end            
            if length(iChan) ~= length(transformPath)
                error('A sparate path must be specified for each channel!')
            end

            for j = 1:length(iChan)

                if exist(transformPath{j},'file')
                    obj.funParams_.TransformFilePaths(iChan(j)) = transformPath(j);
                else
                   error(['The transform file name specified for channel ' ...
                       num2str(iChan(j)) ' is not valid!!']) 
                end
            end
            
        end
        
        function transforms = getTransformation(obj,iChan)                                                                        
            
            %Loads and checks specified transformation(s).            
            if ~ obj.checkChanNum(iChan)
                error('Invalid channel index!')
            end
            
            nChan = length(iChan);
            transforms = cell(1,nChan);
            
            for j = 1:nChan   
                
                tmp = load(obj.funParams_.TransformFilePaths{iChan(j)});
                
                fNames = fieldnames(tmp);
                
                isXform = cellfun(@(x)(istransform(tmp.(x))),fNames);
                
                if ~any(isXform)
                    error(['The transform file specified for channel ' ...
                        num2str(iChan(j)) ...
                        '  does not contain a valid image transformation!']);
                elseif sum(isXform) > 1
                    error(['The transform file specified for channel ' ...
                        num2str(iChan(j)) ...
                        '  contains more than one valid image transformation!']);
                else                
                    transforms{j} = tmp.(fNames{isXform});
                end
                    
                
            end
            
        end        
        
    end
    methods(Static)
        function name =getName()
            name = 'Mask Transformation';
        end
    end
end
    