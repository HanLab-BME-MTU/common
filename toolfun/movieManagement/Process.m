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
        procChanged_ % procChanged_ = true if process is (1) outdated OR 
                     % (2) the process is up to date but process's para 
                     % are changed
                     %
                     % (1) If parent process is changed (true), then the
                     %  children processes will also be set to 'true'
                     % after exception array is returned by sanitycheck   
                     % 
                     % (2) Whether existing process's pamameters are changed 
                     % since last run. (if the process is outdated)
                     % changed - true, unchanged - false

                     
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
        
        % Temp function
        function runProcess(obj) % throws exception
            pause(.5);
            if abs(normrnd(1, 2)) > 100
                error('lccb:runtime:fatal',['T_T Runtime error occurs in step' ...
                    obj.name_ '... zz ZZZ\n\n']);
            end
        end
    end
    methods (Abstract)
        sanityCheck(obj)
        % More abstract classed goes here
        % ... ...
    end
end