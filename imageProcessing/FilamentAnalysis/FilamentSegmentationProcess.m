classdef FilamentSegmentationProcess < ImageProcessingProcess
    % A concrete class for steerable filtering
    %
    % Liya Ding, 06. 2012
    
    methods (Access = public)
        
        function obj = FilamentSegmentationProcess(owner,name,funName,funParams,...
                inImagePaths,outImagePaths, FilamentSegmentationImagePaths)
            
            if nargin == 0
                super_args = {};
            else
                
                super_args{1} = owner;
                super_args{2} = name;
                if nargin > 2
                    super_args{3} = funName;
                end
                if nargin > 3
                    super_args{4} = funParams;
                end
                
                super_args{3} = @filament_segmentation;
                
                if nargin < 3 || isempty(funParams)
                    if nargin <2, outputDir = owner.outputDirectory_; end
                    funParams=FilamentSegmentationProcess.getDefaultParams(owner,outputDir);
                end
                
                if nargin > 4
                    super_args{5} = inImagePaths;
                end
                if nargin > 5
                    super_args{6} = outImagePaths;
                end
                
            end
            
            obj = obj@ImageProcessingProcess(super_args{:});
            
            if nargin > 6
                obj.inFilePaths_(2,:) = FilamentSegmentationImagePaths;
            else
                obj.inFilePaths_(2,:) = cell(1,numel(owner.channels_));
            end
            
        end
        
        function setFilamentSegmentationImagePath(obj,iChan,imagePaths)
            if ~obj.checkChanNum(iChan);
                error('lccb:set:fatal','Invalid image channel number for correction image path!\n\n');
            end
            nChan = length(iChan);
            if ~iscell(imagePaths)
                imagePaths = {imagePaths};
            end
            if numel(imagePaths) ~= nChan
                error('lccb:set:fatal','You must specify one image path for each correction image channel!\n\n');
            end
            for j = 1:nChan
                if exist(imagePaths{j},'dir') && numel(imDir(imagePaths{j})) > 0
                    obj.inFilePaths_{2,iChan(j)} = imagePaths{j};
                else
                    error(['The correction image path specified for channel ' num2str(iChan(j)) ' was not a valid image-containing directory!'])
                end
            end
        end
        
        function fileNames = getFilamentSegmentationImageFileNames(obj,iChan)
            if obj.checkChanNum(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.inFilePaths_(2,iChan),'FilamentSegmentation',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'FilamentSegmentation',false)),fileNames,'FilamentSegmentation',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if any(nIm == 0)
                    error('No images in one or more orientation channels!')
                end
            else
                error('lccb:set:fatal','Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
        end
        
        function setOutImagePath(obj,chanNum,imagePath)
            
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
                    obj.outFilePaths_{1,chanNum(j)} = imagePath{j};
                end
            end
            
            
        end
        function h = draw(obj,iChan,varargin)
            
            outputList = obj.getDrawableOutput();
            drawFilamentSegmentationImage = any(strcmpi('FilamentSegmentationImage',varargin));
            
            if drawFilamentSegmentationImage
                % Input check
                ip =inputParser;
                ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
                ip.addParamValue('output',[],@ischar);
                ip.KeepUnmatched = true;
                ip.parse(iChan,varargin{:})
                
                % Load average corrected image
                s = load(obj.outFilePaths_{2,iChan});
                tmpFields=fieldnames(s);
                data=s.(tmpFields{1});
                
                iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
                if ~isempty(outputList(iOutput).formatData),
                    data=outputList(iOutput).formatData(data);
                end
                
                try
                    assert(~isempty(obj.displayMethod_{iOutput,iChan}));
                catch ME
                    obj.displayMethod_{iOutput,iChan}=...
                        outputList(iOutput).defaultDisplayMethod(iChan);
                end
                
                % Delegate to the corresponding method
                tag = [obj.getName '_channel' num2str(iChan) '_output' num2str(iOutput)];
                drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                    2*numel(fieldnames(ip.Unmatched)),1);
                h=obj.displayMethod_{iOutput,iChan}.draw(data,tag,drawArgs{:});
            else
                h=draw@ImageProcessingProcess(obj,iChan,varargin{1},varargin{2:end});
            end
        end
    end

    
    methods (Static)
        function name =getName()
            name = 'Steerable filtering';
        end
        function h = GUI()
            h= @filamentSegmentationProcessGUI;
        end
        
        function output = getDrawableOutput()
            output = ImageProcessingProcess.getDrawableOutput();
            output(2).name='Steerable Filtering images';
            output(2).var='SteerableFilteringImage';
            output(2).formatData=[];
            output(2).type='graph';
            output(2).defaultDisplayMethod=@(x)LineDisplay('Color',[0 0 0],...
                'LineStyle','-','LineWidth',2,...
                'XLabel','Frame Number','YLabel','FilamentSegmentationImage');
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory =  [outputDir  filesep 'FilamentSegmentationImage'];
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.MaskChannelIndex = funParams.ChannelIndex;
            funParams.BatchMode = false;
        end
    end
end