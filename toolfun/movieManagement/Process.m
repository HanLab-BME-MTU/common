classdef Process < handle 
% Defines the abstract class Process from which every user-defined process 
% will inherit.
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
       name_ 
       createTime_ % Put down the time when process object is created
    end
    properties(SetAccess = protected, GetAccess = public)
        owner_
        hasChanged_ 
        notes_
    end
    methods (Access = protected)
        function obj = Process(owner, name)
            % Constructor of class Process
            if nargin > 0
                obj.owner_ = owner; 
                obj.name_ = name;
                obj.createTime_ = clock;
                obj.hasChanged_ = false;
            end
        end
        function setNotes(obj, text)
            obj.notes_ = text;
        end
    end
    methods (Abstract)
        sanityCheck(obj)

        % More abstract classed goes here
        % ... ...
    end
end