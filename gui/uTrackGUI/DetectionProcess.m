classdef DetectionProcess < ImageAnalysisProcess
% A class definition for a generic detection process.
%
% Chuangang Ren, 11/2010
 
properties(SetAccess = protected, GetAccess = public)
    
    channelIndex_ % The index of channel to process
    overwrite_ = 0; % If overwrite the original result MAT file
end

methods(Access = public)
    
    function obj = DetectionProcess(owner, name, funName, channelIndex, funParams )
       
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
           obj.channelIndex_ = channelIndex; 
        end
        
        if nargin > 4
            obj.funParams_ = funParams;
        end
        

    end
    
    function sanityCheck(obj) 
    end
        
    % Set overwrite
    function setOverwrite (obj, i)
       obj.overwrite_ = i;
    end    

    
    function setChannelIndex(obj, index)
        % Set channel index
        if any(index > length(obj.owner_.channels_))
           error ('User-defined: channel index is larger than the number of channels.') 
        end
        obj.channelIndex_ = index;
    end
    
    function hfigure = resultDisplay(obj,fig,procID)
        % Display the output of the process
        
        % Copied and pasted from the old uTrackPackageGUI
        % but there is definitely some optimization to do
        % Check for movie output before loading the GUI
        chan = [];
        for i = 1:length(obj.owner_.channels_)
            if obj.checkChannelOutput(i)
                chan = i;
                break
            end
        end
        
        if isempty(chan)
            warndlg('The current step does not have any output yet.','No Output','modal');
            return
        end
        
        % Make sure detection output is valid
        load(obj.outFilePaths_{chan},'movieInfo');
        firstframe = [];
        for i = 1:length(movieInfo)
            
            if ~isempty(movieInfo(i).amp)
                firstframe = i;
                break
            end
        end
        
        if isempty(firstframe)
            warndlg('The detection result is empty. There is nothing to visualize.','Empty Output','modal');
            return
        end
        
        if isa(obj, 'Process')
            hfigure = detectionVisualGUI('mainFig', fig, procID);
        else
            error('User-defined: the input is not a Process object.')
        end
    end
    
end


end