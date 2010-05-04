classdef Process < handle 
% Defines the abstract class Process from which every user-defined process 
% will inherit.
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
       name_ 
       createTime_  % Put down the time when process object is created
    end
    properties(SetAccess = protected, GetAccess = public)
        owner_
        procChanged_ % Whether existing process's pamameters are changed 
                     % since last run
        success_ % If the process is successfully run in the most recent 
                 % time
        notes_
        
       funName_
       funParams_        
    end
    methods (Access = protected)
        function obj = Process(owner, name)
            % Constructor of class Process
            if nargin > 0
                obj.owner_ = owner;
                obj.name_ = name;
                
                obj.createTime_ = clock;
                obj.procChanged_ = false;
                obj.success_ = false;
            end
        end
    end
    methods(Access = public)
        function setSuccess(obj, is)
            % set process status 'is' true or false
            obj.success_ = is;
        end
        function setNotes(obj, text)
            % Set the note of process
            obj.notes_ = text;
        end
        function setPara(obj, para)
            % Reset process' parameters
            obj.funParams_ = para;
        end
        function setProcChanged(obj, changed)
            % Set the status that if process's parameters have been
            % changed. changed 'true' or 'false'
           obj.procChanged_ =  changed;
        end
        function change = hasChanged(obj)
            % 'hasChanged' returns true if satisfying the following two
            % conditions:
            % 1. Process has been successfully run (obj.success_ = ture)
            % 2. Pamameters are changed (reported by uicontrols in setting
            % panel, and obj.procChanged_ field is 'true')
            if obj.success_ && obj.procChanged_
                change = true;
            else
                change = false;
            end
        end        
        
    end
    methods (Abstract)
        sanityCheck(obj)
        % More abstract classed goes here
        % ... ...
    end
end