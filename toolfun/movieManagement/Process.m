classdef Process < handle 
% Defines the abstract class Process from which every user-defined process 
% will inherit.
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
       name_ 
       owner_
       timeStamp_ % Put down the time when process object is created
    end
    
    methods
        function obj = Process(owner, name)
            % Constructor of class Process
            if nargin > 0
                obj.owner_ = owner; 
                obj.name_ = name;
                obj.timeStamp_ = clock;
            end
        end      
    end
    methods (Abstract)
        % Abstract methods to be implemented by sub-classes
        isValid = sanityCheck(obj)
        % More abstract classed goes here
        % ... ...
    end
end