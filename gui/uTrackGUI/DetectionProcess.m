classdef DetectionProcess < ImageAnalysisProcess
    % A class definition for a generic detection process.
    %
    % Chuangang Ren, 11/2010
    % Sebastien Besson (last modified Dec 2011)
    
    methods(Access = public)
        
        function obj = DetectionProcess(owner, name, funName, funParams )
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;
            end
            if nargin > 3
                obj.funParams_ = funParams;
            end
        end
        
        function status = checkChannelOutput(obj,iChan)
            
            %Checks if the selected channels have valid output files
            nChan = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan), iChan = 1:nChan; end
            
            status=  ismember(iChan,1:nChan) & ....
                arrayfun(@(x) exist(obj.outFilePaths_{1,x},'file'),iChan);
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'movieInfo'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            s = load(obj.outFilePaths_{iChan},output{:});
           
            if numel(ip.Results.iFrame)>1,
                varargout{1}=s.(output{1});
            else
                varargout{1}=s.(output{1})(iFrame);
            end
        end
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Objects';
            output(1).var='movieInfo';
            output(1).formatData=@formatDetectionOutput;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none','Color',colors(x,:));
        end 
        
        function hfigure = resultDisplay(obj,fig,procID)
            % Display the output of the process
              
            % Check for movie output before loading the GUI
            iChan = find(obj.checkChannelOutput,1);         
            if isempty(iChan)
                warndlg('The current step does not have any output yet.','No Output','modal');
                return
            end
            
            % Make sure detection output is valid
            movieInfo=obj.loadChannelOutput(iChan,'output','movieInfo');
            firstframe=find(arrayfun(@(x) ~isempty(x.amp),movieInfo),1);
            if isempty(firstframe)
                warndlg('The detection result is empty. There is nothing to visualize.','Empty Output','modal');
                return
            end
            
            hfigure = detectionVisualGUI('mainFig', fig, procID);
        end
        
        
    end
    methods(Static)
        function name = getName()
            name = 'Detection';
        end
        function h = GUI()
            h= @detectionProcessGUI;
        end
        function procClasses = getConcreteClasses()
            procClasses = ...
                {'SubResolutionProcess';
                'NucleiDetectionProcess'};
        end
    end

end

function y =formatDetectionOutput(x)
if isempty(x.xCoord)
    y=NaN(1,2);
else
    y = horzcat(x.yCoord(:,1),x.xCoord(:,1));
end
end