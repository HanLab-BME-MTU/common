
classdef Process < hgsetget
    % Defines the abstract class Process from which every user-defined process
    % will inherit.
    %
    
    properties (SetAccess = private, GetAccess = public)
        name_           % Process name
        owner_          % Movie data object owning the process
        createTime_     % Time process was created  
        startTime_      % Time process was last started
        finishTime_     % Time process was last run
    end
    
    properties  (SetAccess = protected)
        % Success/Uptodate flags
        procChanged_   % Whether process parameters have been changed     
        success_       % If the process has been successfully run
        % If the parameter of parent process is changed
        % updated_ - false, not changed updated_ - true
        updated_ 
        
        funName_        % Function running the process
        funParams_      % Parameters for running the process
        visualParams_   % Visualization parameters    
        
        inFilePaths_    % Path to the process input
        outFilePaths_   % Path to the process output

    end
    properties        
        notes_          % Process notes
    end
    properties (Transient=true)
        displayMethod_  % Cell array of display methods
    end
    methods (Access = protected)
        function obj = Process(owner, name)
            % Constructor of class Process
            if nargin > 0
                %Make sure the owner is a MovieData object
                if isa(owner,'MovieObject')
                    obj.owner_ = owner;
                else
                    error('lccb:Process:Constructor','The owner of a process must always be a movie object!')
                end
                
                if nargin > 1
                    obj.name_ = name;
                end
                obj.createTime_ = clock;
                obj.procChanged_ = false;
                obj.success_ = false;
                obj.updated_ = true;
            end
        end
    end
    
    methods
        
        function setPara(obj, para)
            % Reset process' parameters
            if ~isequal(obj.funParams_,para)
                obj.funParams_ = para;
                obj.procChanged_=true;
            end
        end
        
        function setVisualParams(obj, para)
            obj.visualParams_ = para;
        end
        
        function setUpdated(obj, is)
            % Set update status of the current process
            % updated - true; outdated - false
            obj.updated_ = is;
        end
        
        function setDateTime(obj)
            %The process has been re-run, update the time.
            obj.finishTime_ = clock;
        end
        
        function status = checkChanNum(obj,iChan)
            %Checks that the selected channel numbers are valid for this
            %movieData. Specific processes can override this to add more checks.
            if nargin < 2 || isempty(iChan)
                error('You must specify a channel number!')
            end
            
            status = ismember(iChan,1:numel(obj.owner_.channels_));
        end
        
        function status = checkFrameNum(obj,iFrame)
            if nargin < 2 || isempty(iFrame)
                error('You must specify a frame number!')
            end
            
            status = ismember(iFrame,1:obj.owner_.nFrames_);
        end
   
        function sanityCheck(obj)
            crtParams=obj.funParams_;
            defaultParams = obj.getDefaultParams(obj.owner_);
            
            crtFields = fieldnames(crtParams);
            defaultFields = fieldnames(defaultParams);
            status = ~ismember(defaultFields,crtFields);
            for i=find(status)'
                 crtParams.(defaultFields{i})=defaultParams.(defaultFields{i});
            end
            obj.setPara(crtParams);

        end
        
        function run(obj,varargin)
            % Run the process!
            obj.success_=false;
            obj.startTime_ = clock;
            try
                obj.funName_(obj.owner_,varargin{:});
            catch runException
                rethrow(runException)
            end
            obj.success_=true;
            obj.updated_=true;
            obj.procChanged_=false;
            obj.finishTime_ = clock;
            obj.owner_.save;
        end
        
        
        function setInFilePaths(obj,paths)
            %  Set input file paths
            obj.inFilePaths_=paths;
        end
        
        function setOutFilePaths(obj,paths)
            % Set output file paths
            obj.outFilePaths_ = paths;
        end
        
        function time = getProcessingTime(obj)
            %The process has been re-run, update the time.
            time=sec2struct(24*3600*(datenum(obj.finishTime_)-datenum(obj.startTime_)));
        end
        
        function [packageID procID] = getPackage(obj)
            % Retrieve package to which the process is associated
            isOwner=@(x)cellfun(@(y) isequal(y,obj),x.processes_);
            validPackage = cellfun(@(x) any(isOwner(x)),obj.owner_.packages_);
            packageID = find(validPackage);
            procID= cellfun(@(x) find(isOwner(x)),obj.owner_.packages_(validPackage));
        end
        
        function relocate(obj,oldRootDir,newRootDir)
            % Relocate all paths in various fields of process
            %
            % Sebastien Besson, 5/2011
            
            relocateFields ={'inFilePaths_','outFilePaths_',...
                'funParams_','visualParams_'};
            for i=1:numel(relocateFields)
                obj.(relocateFields{i}) = relocatePath(obj.(relocateFields{i}),...
                    oldRootDir,newRootDir);
            end
            
        end
        
        function hfigure = resultDisplay(obj)
            hfigure = movieDataViewer(obj.owner_, ...
                find(cellfun(@(x)isequal(x,obj),obj.owner_.processes_)));    
        end
        
        function h=draw(obj,iChan,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iChan',@isnumeric);
            ip.addRequired('iFrame',@isnumeric);
            ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj,iChan,iFrame,varargin{:})
			
            data=obj.loadChannelOutput(iChan,iFrame,'output',ip.Results.output);
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
        end
        
    end

    methods (Static,Abstract)
        getDefaultParams
        getName
    end
end