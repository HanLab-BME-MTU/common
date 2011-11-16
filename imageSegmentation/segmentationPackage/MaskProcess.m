classdef MaskProcess < Process
    % An abstract class for mask process
    %
    % Sebastien Besson Oct 2011
    
    methods (Access = protected)
        function obj = MaskProcess(owner,name,funName, funParams,outFilePaths)
            % Constructor of class MaskProcess
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
                if ~isempty(outFilePaths) && numel(outFilePaths) ...
                        ~= numel(owner.channels_) || ~iscell(outFilePaths)
                    error('lccb:set:fatal','User-defined: Mask paths must be a cell-array of the same size as the number of image channels!\n\n');
                end
                obj.outFilePaths_ = outFilePaths;
                
            else
                obj.outFilePaths_ = cell(1,numel(owner.channels_));
            end
        end
    end
    methods

        %Checks if a particular channel has masks
        function maskStatus = checkChannelOutput(obj,iChan)
            
            nChanTot = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan)
                iChan = 1:nChanTot; %Default is to check all channels
            end
            nChan = numel(iChan);
            maskStatus = false(1,nChan);
            if all(obj.checkChanNum(iChan))
                for j = 1:nChan
                    %Check the directory and number of masks in each
                    %channel.
                    if exist(obj.outFilePaths_{iChan(j)},'dir') && ...
                            length(imDir(obj.outFilePaths_{iChan(j)})) == obj.owner_.nFrames_;
                        maskStatus(j) = true;
                    end
                end
            end
            
        end
        function setOutMaskPath(obj,chanNum,maskPath)
            if obj.checkChanNum(chanNum)
                obj.outFilePaths_{1,chanNum} = maskPath;
            else
                error('lccb:set:fatal','Invalid mask channel number for mask path!\n\n');
            end
        end
        function fileNames = getOutMaskFileNames(obj,iChan)
            if obj.checkChanNum(iChan)
                fileNames = cellfun(@(x)(imDir(x)),obj.outFilePaths_(iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)
                    error('Incorrect number of masks found in one or more channels!')
                end
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
            
            
        end
        
        function mask = loadChannelOutput(obj,iChan,iFrame,varargin)      
            % Input check
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addRequired('iFrame',@(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',[],@ischar);            
            ip.parse(obj,iChan,iFrame,varargin{:})

            
            % Data loading
            maskNames = obj.getOutMaskFileNames(iChan);
            mask =imread([obj.outFilePaths_{iChan} filesep maskNames{1}{iFrame}]);
%             mask=cell(size(iChan));
%             for i=iChan
%                 maskNames = obj.getOutMaskFileNames(i);
%                 mask{i} = arrayfun(@(j) imread([obj.outFilePaths_{i} filesep...
%                     maskNames{1}{j}]),iFrame,'Unif',0);
%             end
        end
        function output = getDrawableOutput(obj)
            output(1).name='Masks';
            output(1).var='';
            output(1).formatData=@getMaskBoundaries;
            output(1).type='overlay';
            colors = hsv(numel(obj.owner_.channels_));
            output(1).defaultDisplayMethod=@(x) LineDisplay('Color',colors(x,:));
        end
        
    end
    methods(Static)
        function procClasses = getConcreteClasses()
            procClasses = ...
                {'ThresholdProcess';
                'MSSSegmentationProcess'};
        end
        
    end
end

function boundaries = getMaskBoundaries(mask)

b=bwboundaries(mask);
b2 =cellfun(@(x) vertcat(x,[NaN NaN]),b,'Unif',false);
boundaries =vertcat(b2{:});

end
