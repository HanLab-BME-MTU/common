classdef Process < hgsetget
    % Defines the abstract class Process from which every user-defined process
    % will inherit.
    %
    
    properties (SetAccess = private, GetAccess = public)
        name_           % Process name
        owner_          % Movie data object owning the process
        dateTime_       % Time process was last run
    end
    
    properties  (SetAccess = protected)
        % Whether process parameter has been changed
        %
        % changed - true, unchanged - false
        procChanged_    
        
        success_ % If the process has been successfully run
        
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
    methods (Access = protected)
        function obj = Process(owner, name)
            % Constructor of class Process
            if nargin > 0
                %Make sure the owner is a MovieData object
                if isa(owner,'MovieData')
                    obj.owner_ = owner;
                else
                    error('lccb:Process:Constructor','The owner of a process must always be a MovieData object!')
                end
                
                if nargin > 1
                    obj.name_ = name;
                end
                obj.dateTime_ = clock;
                obj.procChanged_ = false;
                obj.success_ = false;
                obj.updated_ = true;
            end
        end
    end
    
    methods
        
        function setPara(obj, para)
            % Reset process' parameters
            obj.funParams_ = para;
        end
        
        function setVisualParams(obj, para)
            obj.visualParams_ = para;
        end
        
        function setProcChanged(obj, changed)
            % Set the status that if process's parameters have been
            % changed. changed 'true' or 'false'
            obj.procChanged_ =  changed;
        end
        
        function setUpdated(obj, is)
            % Set update status of the current process
            % updated - true; outdated - false
            obj.updated_ = is;
        end
        
        function setSuccess(obj, is)
            % Set the success status of process
            obj.success_ = is;
        end
        
        
        function setDateTime(obj)
            %The process has been re-run, update the time.
            obj.dateTime_ = clock;
        end
        
        function OK = checkChanNum(obj,iChan)
            %Checks that the selected channel numbers are valid for this
            %movieData. Specific processes can override this to add more checks.
            if nargin < 2 || isempty(iChan)
                error('You must specify a channel number!')
            end
            
            OK = arrayfun(@(x)(x <= numel(obj.owner_.channels_) && ...
                x > 0 && isequal(round(x),x)),iChan);
        end
        
        function OK = checkFrameNum(obj,iFrame)
            if nargin < 2 || isempty(iFrame)
                error('You must specify a frame number!')
            end
            
            OK = arrayfun(@(x)(x <= obj.owner_.nFrames_ && ...
                x > 0 && isequal(round(x),x)),iFrame);
        end
        
        function run(obj)
            % Run the process!
            obj.funName_(obj.owner_ );
        end
        
        function setInFilePaths(obj,paths)
            %  Set input file paths
            obj.inFilePaths_=paths;
        end
        
        
        function setOutFilePaths(obj,paths)
            % Set output file paths
            obj.outFilePaths_ = paths;
        end
    end
    methods (Abstract)
        sanityCheck(obj)
        % More abstract classed goes here
        % ... ...
    end
end