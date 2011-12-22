classdef  MovieObject < hgsetget
    % Interface defining movie analyis tools
    
    properties (SetAccess = protected)
        createTime_             % Time movie object is created
        processes_ = {};        % Process object cell array
        packages_ = {};         % Package object cell array
    end
    
    properties
        outputDirectory_  ='';      % Default output directory for all processes
        notes_                  % User's notes
    end
    
    methods
        %% Set/get methods
        function set.outputDirectory_(obj, path)
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            path = regexprep(path,endingFilesepToken,'');
            obj.checkPropertyValue('outputDirectory_',path);
            obj.outputDirectory_=path;
        end
        
        function checkPropertyValue(obj,property, value)
            % Check if a property/value pair can be set up
            
            % Test if the property is unchanged
            if isequal(obj.(property),value), return; end
            
            if  ~isempty(obj.(property)),
                % Test if the property is writable
                if ~obj.checkProperty(property)
                    propertyName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
                    error(['The ' propertyName ' has been set previously and cannot be changed!']);
                end
            end
            
            % Test if the value is valid
            if ~obj.checkValue(property,value);
                propertyName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
                error(['The supplied ' propertyName ' is invalid!']);
            end
            
        end
        
        function set.notes_(obj, value)
            obj.checkPropertyValue('notes_',path);
            obj.notes_=value;
        end
        
        
        %% Functions to manipulate process object array
        function addProcess(obj, newprocess)
            obj.processes_ = horzcat(obj.processes_, {newprocess});
        end
        
        function deleteProcess(obj, process)
            % Delete and clear given process object in movie data's process array
            %
            % Input:
            %        process - Process object or index of process to be
            %                  deleted in movie data's process list
            
            % Check input
            if isa(process, 'Process')
                pid = obj.getProcessIndex(process,1,Inf,false);
                if isempty(pid)
                    error('User-defined: The given process is not in current movie processes list.')
                elseif length(pid) ~=1
                    error('User-defined: More than one process of this type exists in movie processes list.')
                end
                
                % Make sure process is a integer index
            elseif ismember(process,1:numel(obj.processes_))
                pid = process;
            else
                error('Please provide a Process object or a valid process index of movie data processes list.')
            end
            
            % Unassociate process in corresponding packages
            if ~isempty(obj.processes_{pid})
                [packageID procID] = obj.processes_{pid}.getPackage;
                if ~isempty(packageID)
                    for i=1:numel(packageID)
                        obj.packages_{packageID(i)}.setProcess(procID(i),[]);
                    end
                end
            end
            
            % Delete and clear the process object
            delete(obj.processes_{pid})
            obj.processes_(pid) = [ ];
        end
        
        function replaceProcess(obj, pid, newprocess)
            % Input check
            ip=inputParser;
            ip.addRequired('obj');
            ip.addRequired('pid',@(x) isscalar(x) && ismember(x,1:numel(obj.processes_)) || isa(x,'Process'));
            ip.addRequired('newprocess',@(x) isa(x,'Process'));
            ip.parse(obj, pid, newprocess);
            
            % Retrieve process index if input is of process type
            if isa(pid, 'Process')
                pid = find(cellfun(@(x)(isequal(x,pid)),obj.processes_));
                assert(isscalar(pid))
            end
            
            % Check new process is compatible with the parent package
            [packageID procID] = obj.processes_{pid}.getPackage;
            if ~isempty(packageID)
                for i=1:numel(packageID)
                    checkNewProcClass = isa(newprocess,...
                        obj.packages_{packageID(i)}.getProcessClassNames{procID(i)});
                    if ~checkNewProcClass
                        error('Package compatibility prevents process replacement');
                    end
                end
            end
            
            % Delete old process and replace it by the new one
            oldprocess=obj.processes_{pid};
            obj.processes_{pid} = newprocess;
            delete(oldprocess);
            if ~isempty(packageID),
                for i=1:numel(packageID)
                    obj.packages_{packageID(i)}.setProcess(procID(i),newprocess);
                end
            end
        end
        
        function iProc = getProcessIndex(obj,procName,varargin)
            % Find the index of a process or processes with given class name
            %
            % SYNOPSIS      iProc=obj.getProcessIndex(procName)
            %               iProc=obj.getProcessIndex(procName,nDesired,askUser)
            
            % Input check
            ip = inputParser;
            ip.addRequired('procName',@ischar);
            ip.addOptional('nDesired',1,@isscalar);
            ip.addOptional('askUser',true,@isscalar);
            ip.parse(procName,varargin{:});
            nDesired = ip.Results.nDesired;
            askUser = ip.Results.askUser;
            
            % Read process of given type
            iProc = find(cellfun(@(x)(isa(x,procName)),obj.processes_));
            nProc = numel(iProc);
            
            %If there are only nDesired or less processes found, return
            if nProc <= nDesired, return; end
            
            % If more than nDesired processes
            if askUser
                isMultiple = nDesired > 1;
                procNames = cellfun(@(x)(x.getName),...
                    obj.processes_(iProc),'UniformOutput',false);
                iSelected = listdlg('ListString',procNames,...
                    'SelectionMode',isMultiple,...
                    'ListSize',[400,400],...
                    'PromptString',['Select the desired ' procName ':']);
                iProc = iProc(iSelected);
                if isempty(iProc)
                    error('You must select a process to continue!');
                end
            else
                warning('lccb:process',['More than ' num2str(nDesired) ' ' ...
                    procName 'es were found! Returning most recent process(es)!'])
                iProc = iProc(end:-1:(end-nDesired+1));
            end
        end
        
        %% Functions to manipulate package object array
        function addPackage(obj, newpackage)
            % Add a package to the packages_ array
            obj.packages_ = horzcat(obj.packages_ , {newpackage});
        end
        
        function deletePackage(obj, package)
            % Check input
            if isa(package, 'Package')
                pid = find(cellfun(@(x)isequal(x, package), obj.packages_));
                if isempty(pid)
                    error('User-defined: The given package is not in current movie processes list.')
                elseif length(pid) ~=1
                    error('User-defined: More than one process of this type exists in movie processes list.')
                end
                
                % Make sure process is a integer index
            elseif ismember(package,1:numel(obj.packages_))
                pid = package;
            else
                error('Please provide a Package object or a valid package index of movie data processes list.')
            end
            
            % Delete and clear the process object
            delete(obj.packages_{pid})
            obj.packages_(pid) = [ ];
        end
        
        %% Miscellaneous functions
        
        function [askUser,relocateFlag] = sanityCheck(obj, path, filename,askUser)
            % Check if the path and filename stored in the movie object are
            % the same as the ones provided in argument. They can differ if
            % the MAT file has been renamed or moved to another location.
            
            if nargin < 4, askUser = true; end
            if nargin > 1
                % Remove ending file separators from paths
                endingFilesepToken = [regexptranslate('escape',filesep) '$'];
                oldPath = regexprep(obj.getPath(),endingFilesepToken,'');
                newPath = regexprep(path,endingFilesepToken,'');
                if ~strcmp(oldPath, newPath)
                    if askUser
                        objName = regexprep(class(obj),'([A-Z])',' ${lower($1)}');
                        relocateMsg=sprintf(['The' objName ' located in \n%s\n has been relocated to \n%s.\n'...
                            'Should I try to relocate the components of the' objName ' as well?'],oldPath,newPath);
                        confirmRelocate = questdlg(relocateMsg,['Relocation -' objName],'Yes to all','Yes','No','Yes');
                    else
                        confirmRelocate = 'Yes to all';
                    end
                    
                    % Relocate
                    relocateFlag = ~strcmp(confirmRelocate,'No');
                    if relocateFlag
                        obj.relocate(path);
                        askUser = strcmp(confirmRelocate,'Yes');
                    else
                        obj.setPath(path);
                    end
                end
                obj.setFilename(filename);
            end
        end
        
        function [oldRootDir newRootDir]=relocate(obj,newPath)
            % Relocate all paths of the movie object
            %
            % This function automatically relocates the output directory,
            % processes and package paths assuming the internal
            % architecture of the  project is conserved.
            
            [oldRootDir newRootDir]=getRelocationDirs(obj.getPath,newPath);
            
            obj.outputDirectory_=relocatePath(obj.outputDirectory_,oldRootDir,newRootDir);
            obj.setPath(newPath);
            
            % Relocate paths in processes input/output as well as function
            % and visual parameters
            for i=1:numel(obj.processes_),
                obj.processes_{i}.relocate(oldRootDir,newRootDir);
            end
            
            for i=1:numel(obj.packages_),
                obj.packages_{i}.relocate(oldRootDir,newRootDir);
            end
        end
        
        function flag = save(obj)
            % Save the movie object to disk.
            %
            % This function check for the values of the path and filename.
            % If empty, it launches a dialog box asking where to save the
            % movie object. If a MAT file already exist, copies this MAT
            % file before saving the MovieObject.
            %
            % OUTPUT:
            %    flag - a flag returning the status of the save method
            movieClass =class(obj);
            
            % If no path or fileName, start a dialog box asking where to save the MovieObject
            if isempty(obj.getPath) || isempty(obj.getFilename)
                if ~isempty(obj.getPath),
                    defaultDir=obj.getPath;
                elseif ~isempty(obj.outputDirectory_)
                    defaultDir=obj.outputDirectory_;
                else
                    defaultDir =pwd;
                end
                objName = regexprep(movieClass,'([A-Z])',' ${lower($1)}');
                defaultName = regexprep(movieClass,'(^[A-Z])','${lower($1)}');
                [filename,path] = uiputfile('*.mat',['Find a place to save your' objName],...
                    [defaultDir filesep defaultName '.mat']);
                
                if ~any([filename,path]), flag=0; return; end
                
                % After checking file directory, set directory to movie data
                obj.setPath(path);
                obj.setFilename(filename);
            end
            
            %First, check if there is an existing file. If so, save a
            %backup. Then save the MovieData as obj
            fullPath = [obj.getPath filesep obj.getFilename];
            if exist(fullPath,'file');
                copyfile(fullPath,[fullPath(1:end-3) 'old']);
            end
            
            switch (movieClass)
                case 'MovieData'
                    MD=obj; %#ok<NASGU>
                    save(fullPath,'MD');
                case 'MovieList'
                    ML=obj; %#ok<NASGU>
                    save(fullPath,'ML');
            end
            
            flag=1;
        end
        
        function reset(obj)
            % Reset the movieObject
            obj.processes_={};
            obj.packages_={};
        end
        
    end
    methods (Abstract)
        getPath(obj)
        getFilename(obj)
        setPath(obj,path)
        setFilename(obj,filename)
    end
    methods(Static)
        function status = checkProperty(property)
            % Returns true/false if the non-empty property is writable
            status = false;
            switch property
                case {'notes_'};
                    status = true;
                case {'outputDirectory_'};
                    stack = dbstack;
                    if any(cellfun(@(x)strcmp(x,'MovieObject.relocate'),{stack.name})),
                        status  = true;
                    end
            end
        end
        
        
        
        function obj = load(moviepath,askUser)
            % Load a movie object stored in a mat file
            
            % Check the path is a valid Mat file
            assert(exist(moviepath, 'file')==2,'lccb:movieObject:load', 'File does not exist.');
            
            if strcmpi(moviepath(end-3:end),'.mat')
                try
                    vars = whos('-file',moviepath);
                catch ME
                    ME2 = MException('lccb:movieObject:load', 'Fail to open file. Make sure it is a MAT file.');
                    ME.addCause(ME2);
                    throw(ME);
                end
                
                % Check if a movie object is part of the v
                isMovie = cellfun(@(x) any(strcmp(superclasses(x),'MovieObject')),{vars.class});
                if ~any(isMovie)
                    error('lccb:movieObject:load', ...
                        'No movie object is found in selected MAT file.');
                end
                if sum(isMovie)>1
                    error('lccb:movieObject:load', ...
                        'Multiple movie objects are found in selected MAT file.');
                end
                
                data = load(moviepath,'-mat',vars(isMovie).name);
                obj= data.(vars(isMovie).name);
                
                if nargin<2, askUser=true; end
                [moviePath,movieName,movieExt]=fileparts(moviepath);
                obj.sanityCheck(moviePath,[movieName movieExt],askUser);
            else
                obj=bfImport(moviepath);
            end
        end
        
    end
end