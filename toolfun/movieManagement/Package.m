classdef Package < handle 
% Defines the abstract class Package from which every user-defined packages 
% will inherit. This class cannot be instantiated.
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % defining accessor functions
    %
    % Objects of sub-class of Package cannot change variable values since 
    % 'SetAccess' attribute is set to private
        name_  % the name of instantiated package
        createTime_ % The time when object is created. Same output as 
                    % function 'clock' 
    end
    properties(SetAccess = protected, GetAccess = public)
        owner_ % The MovieData object this package belongs to
        processes_ % Cell array containing all processes who will be used in this package
        processClassNames_ % Must have accurate process class name
                           % List of processes required by the package, 
                           % Cell array - same order and number of elements
                           % as processes in dependency matrix
        depMatrix_ % Processes' dependency matrix
        notes_ % The notes users put down
    end
    
    methods (Access = protected)
        function obj = Package(owner, name, depMatrix, processClassNames)
            % Constructor of class Package
            
            if nargin > 0
                obj.name_ = name;
                obj.createTime_ = clock;
                obj.owner_ = owner; 
                obj.depMatrix_ = depMatrix;
                
                obj.processClassNames_ = processClassNames;
                obj.processes_ = cell(1,length(processClassNames));
               
            end
        end
        function processExceptions = checkProcesses(obj, full) 
            % checkProcesses is called by package's sanitycheck. It returns
            % a cell array of exceptions.
            % The following things will be checked in this function
            %   1. The process itself has a problem
            %   2. The parameters in the process setting panel have changed
            %   3. The process that current process depends on has a
            %      problem
            %
            % processExceptions is a cell array with same length of
            % processClassNames_. It takes down all the exceptions found in
            % sanity check. Exceptions of i th process will be saved in
            % processExceptions{i}
            if nargin < 2
                full = true;
            end

            nProcesses = length(obj.processClassNames_);
            processExceptions = cell(1,nProcesses);
            processVisited = false(1,nProcesses);
            % 1. Check if the process itself has a problem
            for i = 1:nProcesses
                if isempty(obj.processes_{i})
                    continue;
                else
                    try
                        obj.processes_{i}.sanityCheck;
                    catch ME
                        % Add process exception to the ith process
                        processExceptions{i} = horzcat(processExceptions{i}, ME);
                    end
                end
            end
            % 2. Check if the para in the process setting panel
            %    have changed
            for i = 1:nProcesses
                if isempty(obj.processes_{i})
                    continue;
                else            
                    if obj.processes_{i}.hasChanged()
                        ME = MException('lccb:???',['Parameters of this process have been',...
                            ' changed since previous run']);
                        % Add para exception to the ith process
                        processExceptions{i} = horzcat(processExceptions{i}, ME);                           
                    end
                end 
            end
            % 3. Check if the processes that current process depends on
            %    have problems
            
            for i = 1:nProcesses
                if isempty(obj.processes_{i})
                    continue;
                elseif ~processVisited(i)
                   [processExceptions, processVisited]= ...
                            obj.dfs_(i, processExceptions, processVisited);
                end
            end
            % Throw exception array
            throw(processExceptions);
        end
    end
    methods (Access = private)
        function [processExceptions, processVisited] = dfs_(obj, ...
                iProcess,processExceptions,processVisited)
            processVisited(iProcess) = true;
            parentIndex = find(obj.depMatrix_(iProcess,:));
            if isempty(parentIndex)
                return;
            else
                for j = parentIndex
                    if ~isempty(obj.processes_{j}) && ~processVisited(j)
                         [processExceptions, processVisited] = ...
                        obj.dfs_(j, processExceptions,processVisited);
                    end
                    % If j th process has an exception, add an exception to
                    % the exception list of i th process. Since i th
                    % process depends on the j th process
                    if ~isempty(processExceptions{j}) || ...
                                                isempty(obj.processes_{j})
                        ME = MException('lccb:????', ...
                            ['Process ',obj.processClassNames_{iProcess},...
                            'depends on process ',obj.processClassNames_{j},
                            ', which has a problem']);
                        % Add dependency exception to the ith process
                        processExceptions{iProcess} = ...
                                horzcat(processExceptions{iProcess}, ME); 
                    end
                end
            end
        end
    end
    methods (Access = public)
        function setProcess (obj, i, newProcess)
            % set the i th process of obj.processes_ to newprocess
            assert(i<=length(obj.processClassNames_),'UserDefined Error: i exceeds obj.processes length');
            assert(isa(newProcess, 'Process'), ...
                'UserDefined Error: input is not ''Process'' object.');
            obj.processes_{i} = newProcess;
        end
        function setNotes (obj, text)
            obj.notes_ = text;
        end
    end
    
    methods (Abstract)
        sanityCheck(obj)

        % More abstract methods go here
    end
end