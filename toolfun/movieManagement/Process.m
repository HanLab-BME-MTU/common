classdef Process < handle 
% Defines the abstract class Process from which every user-defined process 
% will inherit.
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions              
       dateTime_  % Time process was last run
    end
    properties(SetAccess = protected, GetAccess = public)
        name_  %Process name
        owner_
        procChanged_ % Whether process parameter has been changed 
                     % 
                     % changed - true, unchanged - false
        
        
                     
        success_ % If the process has been successfully run 
        updated_ % If the parameter of parent process is changed
                 % updated_ - false, not changed updated_ - true
        notes_
        
        funName_
        funParams_        
    end
    methods (Access = protected)
        function obj = Process(owner, name)
            % Constructor of class Process
            if nargin > 0
                obj.owner_ = owner;
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
        function setUpdated(obj, is)
            % Set update status of the current process
            % updated - true; outdated - false
            obj.updated_ = is;
        end
        function setDateTime(obj)
            %The process has been re-run, update the time.
            obj.dateTime_ = clock;
        end
        
        % Temp function
%         function runProcess(obj) % throws exception
%             pause(.5);
%             if abs(normrnd(1, 2)) > 10
%                 error('lccb:runtime:fatal',['T_T Runtime error occurs in step' ...
%                     obj.name_ '... zz ZZZ\n\n']);
%             end
%
%            obj.funName_(obj.owner_)
%         end
    end
    methods (Abstract)
        sanityCheck(obj)
        % More abstract classed goes here
        % ... ...
    end
end