classdef Package < hgsetget 
    % Defines the abstract class Package from which every user-defined packages
    % will inherit. This class cannot be instantiated.
    
    properties(SetAccess = immutable)
        createTime_ % The time when object is created.
    end
    
    properties (SetAccess = private)
    % Objects of sub-class of Package cannot change variable values since 
    % 'SetAccess' attribute is set to private
        name_  % the name of instantiated package
    end
 
    properties(SetAccess = protected)
        owner_ % The MovieData object this package belongs to
        processes_ % Cell array containing all processes who will be used in this package
        processClassHandles_ % 
        processClassNames_ % Must have accurate process class name
                           % List of processes required by the package, 
                           % Cell array - same order and number of elements
                           % as processes in dependency matrix
        depMatrix_ % Processes' dependency matrix
        tools_ % Array of external tools
        
    end

    properties
        notes_ % The notes users put down
        outputDirectory_ %The parent directory where results will be stored.
                         %Individual processes will save their results to
                         %sub-directories of this directory.
    end

    methods
        function set.outputDirectory_(obj,value)
            
            if isequal(obj.outputDirectory_,value), return; end
            if ~isempty(obj.outputDirectory_), 
                stack = dbstack;
                if strcmp(stack(3).name,'MovieData.relocate'),
                    error(['This channel''s ' propertyName ' has been set previously and cannot be changed!']);
                end
            end
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            value = regexprep(value,endingFilesepToken,'');
            obj.outputDirectory_=value;
        end
    end
    methods (Access = protected)
        function obj = Package(owner, name, depMatrix, processClassNames, ...
                outputDirectory,varargin)
            % Constructor of class Package
            
            if nargin > 0
                obj.name_ = name;
                obj.owner_ = owner; 
                obj.depMatrix_ = depMatrix;
                obj.processClassNames_ = processClassNames;
                obj.outputDirectory_ = outputDirectory;
                
                nVarargin = numel(varargin);
                if nVarargin > 1 && mod(nVarargin,2)==0
                    for i=1 : 2 : nVarargin-1
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
            
                obj.processes_ = cell(1,length(processClassNames));
                obj.createTime_ = clock;
            end
        end
        
        function [processExceptions, processVisited] = dfs_(obj, ...
                i,processExceptions,processVisited)
            processVisited(i) = true;
            parentIndex = find(obj.depMatrix_(i,:));
            if isempty(parentIndex), return;  end
            for j = parentIndex
                if ~isempty(obj.processes_{j}) && ~processVisited(j)
                    
                    [processExceptions, processVisited] = ...
                        obj.dfs_(j, processExceptions,processVisited);
                end
                % If j th process has an exception, add an exception to
                % the exception list of i th process. Since i th
                % process depends on the j th process
                % Exception is created when satisfy:
                % 1. Process is successfully run in the most recent time
                % 2. Parent process has error OR parent process does
                %    not exist
                if obj.processes_{i}.success_ && ...
                        ( ~isempty(processExceptions{j}) || isempty(obj.processes_{j}) )
                    
                    % Set process's updated=false
                    obj.processes_{i}.setUpdated (false);
                    
                    % Create a dependency error exception
                    ME = MException('lccb:depe:warn', ...
                        ['The current step is out of date because step ',num2str(j),'. ',obj.processes_{j}.name_,', which the current step depends on, is out of date.'...
                        'Please run again to update your result.']);
                    % Add dependency exception to the ith process
                    processExceptions{i} = horzcat(processExceptions{i}, ME);
                end
            end
            
        end
    end
    methods (Access = public)
        
        function processExceptions = sanityCheck(obj, varargin)
            % sanityCheck is called by package's sanitycheck. It returns
            % a cell array of exceptions. Keep in mind, make sure all process
            % objects of processes checked in the GUI exist before running
            % package sanitycheck. Otherwise, it will cause a runtime error
            % which is not caused by algorithm itself.
            %
            % The following steps will be checked in this function
            %   I. The process itself has a problem
            %   II. The parameters in the process setting panel have changed
            %   III. The process that current process depends on has a
            %      problem
            %
            % OUTPUT:
            %   processExceptions - a cell array with same length of
            % processClassNames_. It collects all the exceptions found in
            % sanity check. Exceptions of i th process will be saved in
            % processExceptions{i}
            %
            % INPUT:
            %   obj - package object
            %   full - true   check 1,2,3 steps
            %          false  check 2,3 steps
            %   procID - A. Numeric array: id of processes for sanitycheck
            %            B. String 'all': all processes will do
            %                                      sanity check
            %
            
            nProcesses = length(obj.processClassNames_);
            processExceptions = cell(1,nProcesses);
            processVisited = false(1,nProcesses);
            
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj');
            ip.addOptional('full',true, @(x) islogical(x));
            ip.addOptional('procID',1:nProcesses,@(x) (isvector(x) && ~any(x>nProcesses)) || strcmp(x,'all'));
            ip.parse(obj,varargin{:});
            
            full = ip.Results.full;
            procID = ip.Results.procID;
            if strcmp(procID,'all'), procID = 1:nProcesses;end
            
            validProc = procID(~cellfun(@isempty,obj.processes_(procID)));
            if full 
                % I: Check if the process itself has a problem
                %
                % 1. Process sanity check
                % 2. Input directory  
                for i = validProc
                    try
                        obj.processes_{i}.sanityCheck;
                    catch ME
                        % Add process exception to the ith process
                        processExceptions{i} = horzcat(processExceptions{i}, ME);
                    end
                end
            end
            
            % II: Determine the parameters are changed if satisfying the
            % following two conditions:
            % A. Process has been successfully run (obj.success_ = true)
            % B. Pamameters are changed (reported by uicontrols in setting
            % panel, and obj.procChanged_ field is 'true')
            changedProcesses = validProc(cellfun(@(x) x.success_ && x.procChanged_,obj.processes_(validProc)));
            for i = changedProcesses                    
                % Set process's updated=false
                obj.processes_{i}.setUpdated(false);
                % Create an dependency error exception
                ME = MException('lccb:paraChanged:warn',...
                    'The current step is out of date because the channels or parameters have been changed.');
                % Add para exception to the ith process
                processExceptions{i} = horzcat(processExceptions{i}, ME);
            end
            
            % III: Check if the processes that current process depends
            % on have problems
            for i = validProc
                if ~processVisited(i)
                    [processExceptions, processVisited]= ...
                        obj.dfs_(i, processExceptions, processVisited);
                end
            end
        end
        
        function setDepMatrix(obj,row,col,value)
            % row and col could be array
            obj.depMatrix_(row, col) = value;
        end
            
        function setProcess(obj, i, newProcess)
            % set the i th process of obj.processes_ to newprocess
            % If newProcess = [ ], clear the process in package process
            % list
            assert(i<=length(obj.processClassNames_),...
                'UserDefined Error: i exceeds obj.processes length');
            if isa(newProcess,obj.processClassNames_{i}) ||...
                    isempty(newProcess)      
                obj.processes_{i} = newProcess;
            else
                error('User-defined: input should be Process object or empty.')
            end
        end
        
    end
        

    methods(Static,Abstract)
        start
        getDependencyMatrix
        getOptionalProcessId
    end
    
end