classdef BackgroundSubtractionProcess < ImageCorrectionProcess
    
    %A class for performing background subtraction on images using background masks.
    %
    %Hunter Elliott, 5/2010
    %
    
    methods (Access = public)
        
        function obj = BackgroundSubtractionProcess(owner,outputDir,funParams,backgroundMaskPaths,...
                inImagePaths,outImagePaths)
            if nargin == 0
                super_args = {};
            else
                nChan = numel(owner.channels_);
                
                super_args{1} = owner;
                super_args{2} = BackgroundSubtractionProcess.getName;
                super_args{3} = @backgroundSubtractMovie;
                
                if nargin < 3 || isempty(funParams)
                    
                    %----Defaults----%
                    funParams.OutputDirectory = ...
                        [outputDir  filesep 'background_subtracted_images'];
                    funParams.ChannelIndex = 1:nChan;
                    funParams.MaskChannelIndex = funParams.ChannelIndex;
                    funParams.BatchMode = false;
                    
                    
                end
                
                super_args{4} = funParams;
                
                if nargin > 3
                    %Set the correction image paths to the background mask paths
                    %input.
                    super_args{7} = backgroundMaskPaths;
                end
                
                if nargin > 4
                    super_args{5} = inImagePaths;
                end
                
                if nargin > 5
                    super_args{6} = outImagePaths;
                end
                
            end
            
            obj = obj@ImageCorrectionProcess(super_args{:});
        end
        
        function h = draw(obj,iChan,varargin)
            
            outputList = obj.getDrawableOutput();
            drawBkgValues = any(strcmpi('bkgValues',varargin));
            
            if drawBkgValues
                % Input check
                ip =inputParser;
                ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
                ip.addParamValue('output',[],@ischar);
                ip.KeepUnmatched = true;
                ip.parse(iChan,varargin{:})
                
                % Load average corrected image
                tmp = load(obj.outFilePaths_{2,iChan});
                tmpFields=fieldnames(tmp);
                data(:,1)=tmp.(tmpFields{1});
                data(:,2)=1:size(data,1);
                
                iOutput= 2;
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
            name = 'Background Subtraction';
        end
        function h = GUI()
            h= @backgroundSubtractionProcessGUI;
        end
        
        function output = getDrawableOutput()
            output = ImageProcessingProcess.getDrawableOutput();
            output(2).name='Background values';
            output(2).var='bkgValues';
            output(2).formatData=[];
            output(2).type='graph';
            output(2).defaultDisplayMethod=@(x)LineDisplay('Color',[0 0 0],...
                'XLabel','Frame Number','YLabel','Subtracted Background Value, A.U.');
        end
    end
end
