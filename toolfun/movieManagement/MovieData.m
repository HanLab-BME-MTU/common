classdef  MovieData < hgsetget
% Movie management class used for generic processing
    properties (SetAccess = protected)  
        % Software defined data
    
        nFrames_                % Number of frames
        imSize_                 % Image size 1x2 array[height width]
       
    end
    properties (SetAccess = immutable)
        createTime_             % Time movie data is created
    end
    
    properties
    % User defined data
        
        % ---- Used params ----
        
        channels_ = [];         % Channel object array    
        outputDirectory_        % The default output directory for all processes
        movieDataPath_          % The path where the movie data is saved
        movieDataFileName_      % The name under which the movie data is saved
        notes_                  % User's notes        
        pixelSize_              % Pixel size (nm)
        timeInterval_           % Time interval (s)
        numAperture_            % Numerical Aperture
        camBitdepth_            % Camera Bit-depth
        eventTimes_             % Time of movie events 
        
        % ---- Un-used params ----
        
        magnification_
        binning_
        
    % Progress data
    
        processes_ = {};        % Process object cell array
        packages_ = {};         % Package object cell array

    end

    methods
        %% Constructor
        function obj = MovieData(channels,outputDirectory,varargin)
            % Constructor of the MovieData object
            %
            % INPUT
            %    channels - a Channel object or an array of Channels
            %    outputDirectory - a string containing the output directory
            %    OPTIONAL - a set of options under the property/key format 
            
            % Required input fields
            obj.channels_ = channels;
            obj.outputDirectory_ = outputDirectory;

            % Construct the Channel object
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
            
            obj.createTime_ = clock;
        end

        %% Set/Get methods
        function set.channels_(obj, value)
            obj.checkPropertyValue('channels_',value);
            obj.channels_=value;
        end
        
        function set.outputDirectory_(obj, value)
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            value = regexprep(value,endingFilesepToken,'');
            obj.checkPropertyValue('outputDirectory_',value);
            obj.outputDirectory_=value;
        end
        
        function set.movieDataPath_(obj, value)
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            value = regexprep(value,endingFilesepToken,'');
            obj.checkPropertyValue('movieDataPath_',value);
            obj.movieDataPath_=value;
        end
        
        function set.movieDataFileName_(obj, value)
            obj.checkPropertyValue('movieDataFileName_',value);
            obj.movieDataFileName_=value;
        end     

        function set.notes_ (obj, value)
            obj.checkPropertyValue('notes_',value);
            obj.notes_=value;
        end
        
        function set.pixelSize_ (obj, value)
            obj.checkPropertyValue('pixelSize_',value);
            obj.pixelSize_=value;
        end
        
        function set.timeInterval_ (obj, value)
            obj.checkPropertyValue('timeInterval_',value);
            obj.timeInterval_=value;
        end
        
        function set.numAperture_ (obj, value)
            obj.checkPropertyValue('numAperture_',value);
            obj.numAperture_=value;
        end
        
        function set.camBitdepth_ (obj, value)
            obj.checkPropertyValue('camBitdepth_',value);
            obj.camBitdepth_=value;
        end
        
        function set.magnification_ (obj, value)
            obj.checkPropertyValue('magnification_',value);
            obj.magnification_=value;
        end
        function set.binning_ (obj, value)
            obj.checkPropertyValue('binning_',value);
            obj.binning_=value;
        end

        function checkPropertyValue(obj,property, value)
            % Check if a property/value pair can be set up
            % 
            % Returns an error if either the property is unchangeable or
            % the value is invalid.
            %
            % INPUT:
            %    property - a valid MovieData property name (string)
            %
            %    value - the property value to be checked
            %
            % Sebastien Besson, 4/2011
            
            % Test if the property is writable
            propertyCheck =0;
            if strcmp(property,{'notes_'}), propertyCheck=1;
            elseif isempty(obj.(property)), propertyCheck=1;
            elseif isequal(obj.(property),value), return;
            elseif any(strcmp(property,{'outputDirectory_','movieDataPath_','movieDataFileName_'}))
                % Allow relocation of outputDirectory_, movieDataPath_,'movieDataFileName_'
                stack = dbstack;
                if any(strcmp(stack(3).name,{'MovieData.relocate','MovieData.sanityCheck'})), 
                    propertyCheck=1; 
                end
            end
            
            if ~propertyCheck
                propertyName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
                error(['The ' propertyName ' has been set previously and cannot be changed!']);
            end
            
            % Test if the value is valid
            valueCheck=obj.checkValue(property,value);
            if ~valueCheck
                propertyName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
                error(['The supplied ' propertyName ' is invalid!']);
            end
            
        end
        
        function sanityCheck(obj, movieDataPath, movieDataFileName,askUser)
        % 1. Sanity check (user input, input channels, image files)
        % 2. Assignments to 4 properties:
        %       movieDataPath_
        %       movieDataFileName_
        %       nFrames_
        %       imSize_
           
            % Ask user by default for relocation
            if nargin < 4, askUser = true; end
            
            % Check if the path and filename stored in the movieData are the same
            % as the ones provided in argument. They can differ if the movieData
            % MAT file has been renamed, move or copy to another location.
            if nargin > 1
                
                %Remove ending file separators if any
                endingFilesepToken = [regexptranslate('escape',filesep) '$'];
                path1 = regexprep(obj.movieDataPath_,endingFilesepToken,'');
                path2 = regexprep(movieDataPath,endingFilesepToken,'');
                if  ~strcmp(path1, path2)
                    
                    if askUser
                        relocateMsg=sprintf(['The movie data located in \n%s\n has been relocated to \n%s.\n'...
                            'Should I try to relocate the components of the movie data as well?'],path1,path2);
                        confirmRelocate = questdlg(relocateMsg,'Movie Data','Yes','No','Yes');
                    else
                        confirmRelocate = 'Yes';
                    end
                    
                    if strcmp(confirmRelocate,'Yes')
                        obj.relocate(movieDataPath); 
                    else
                        obj.movieDataPath_=movieDataPath;
                    end
                end
            
                if  ~strcmp(obj.movieDataFileName_, movieDataFileName)
                    obj.movieDataFileName_ = movieDataFileName; 
                end
            
            end
            width = zeros(1, length(obj.channels_));
            height = zeros(1, length(obj.channels_));
            nFrames = zeros(1, length(obj.channels_));
            
            for i = 1: length(obj.channels_)
                [width(i) height(i) nFrames(i)] = obj.channels_(i).sanityCheck(obj);
            end
            
            assert(max(nFrames) == min(nFrames), ...
                'Different number of frames are detected in different channels. Please make sure all channels have same number of frames.')
            
            assert(max(width)==min(width) && max(height)==min(height), ...
                'Image sizes are inconsistent in different channels.\n\n')

		
            % Define imSize_ and nFrames_;     
            if ~isempty(obj.nFrames_)
                assert(obj.nFrames_ == nFrames(1), 'Record shows the number of frames has changed in this movie.')
            else
                obj.nFrames_ = nFrames(1);
            end
            
            if ~isempty(obj.imSize_)
                assert(obj.imSize_(2) == width(1) && obj.imSize_(1) ==height(1), 'Record shows image size has changed in this movie.')
            else
                obj.imSize_ = [height(1) width(1)];
            end

        end
        
        %% Functions to manipulate process object array
        function addProcess(obj, newprocess)
            % Add a process to the processes_ array
            obj.processes_ = horzcat(obj.processes_, {newprocess});
        end
        
        function deleteProcess(obj, process)
           % Delete and clear given process object in movie data's process array
           %
           % Input:
           %        process - Process object or index of process to be
           %                  deleted in movie data's process list

            if isa(process, 'Process')
                
                id = find(cellfun(@(x)isequal(x, process), obj.processes_));
                if isempty(id)
                    error('User-defined: The given process is not in current movie data''s process list.')

                elseif length(id) ~=1
                   error('User-defined: There should be one identical maskrefinement processes exists in movie data''s process list.') 
                end                
                
            % Make sure process is a integer index   
            elseif isnumeric(process) && process == round(process)...
                    && process>0 && process<=length(obj.processes_)
                
                id = process;
            else
                error('Please provide a Process object or a valid process index of movie data''s process list.')
            end
            
            % Delete and clear the process object
            delete(obj.processes_{id})
            obj.processes_(id) = [ ];
            
        end
            
        function iProc = getProcessIndex(obj,procName,nDesired,askUser)
            %Find the index of a process or processes with given class name
            
            if nargin < 2 || isempty(procName)
                error('You must specify a process class name!')
            end            
            if nargin < 3 || isempty(nDesired)
                nDesired = 1;
            end
            if nargin < 4 || isempty(askUser)
                askUser = true;
            end
            if isempty(obj.processes_)
                iProc = [];
                return
            end
            
            iProc = find(arrayfun(@(x)(isa(obj.processes_{x},procName)),1:numel(obj.processes_)));
            nProc = numel(iProc);
            %If there are only nDesired or less processes found, return
            %them.
            if nProc <= nDesired
                return
            elseif nProc > nDesired 
                if askUser
                    isMultiple = nDesired > 1;
                    procNames = cellfun(@(x)(x.name_),...
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
                    warning(['More than ' num2str(nDesired) ' ' ...
                        procName 'es were found! Returning most recent process(es)!'])
                    iProc = iProc(end:-1:(end-nDesired+1));
                end                                
            end
        end
        
        function runProcess(obj, process)
           % Run given process object in movie data's process array
           %
           % Input:
           %        process - Process object or index of process to be
           %                  deleted in movie data's process list

            if isa(process, 'Process')
                
                id = find(cellfun(@(x)isequal(x, process), obj.processes_));
                if isempty(id)
                    error('User-defined: The given process is not in current movie data''s process list.')

                elseif length(id) ~=1
                   error('User-defined: There should be one identical maskrefinement processes exists in movie data''s process list.') 
                end                
                
            % Make sure process is a integer index   
            elseif isnumeric(process) && process == round(process)...
                    && process>0 && process<=length(obj.processes_)
                
                id = process;
            else
                error('Please provide a Process object or a valid process index of movie data''s process list.')
            end
            
            % Run the selected process object
            obj.processes_(id).run();
            
        end

        %% Functions to manipulate package object array
        function addPackage(obj, newpackage)
            % Add a package to the packages_ array
            obj.packages_ = horzcat(obj.packages_ , {newpackage});
        end

        function setPackage(obj, i, package)
            % Set the i th package in packages_ array
            assert( i <= length(obj.packages_), ...
                'Index of package exceeds dimension of packages_ array');
                delete(obj.packages_{i});
                obj.packages_{i} = package;
        end

        
        function relocate(obj,newMovieDataPath)
            % Relocate all paths of the movie data object
            %
            % This function automatically relocates channel and processes 
            % paths assuming the internal architecture of the project is 
            % conserved. 
            % It compares the newMovieDataPath and the stored movieDataPath
            % to guess the new root directory (accouting for OS changes).
            % It then uses the relocatePath function to relocate the
            % channels and output directory paths as well as the internal
            % processes paths (input/ouput, function and visual parameters)
            %
            % Input:
            %
            %   newMovieDataPath - the new location of the movie data
            %   MAT file
            % 
            % Sebastien Besson, 4/2011
            
            %Convert temporarily all path using the local fileseps (for comparison)
            oldMovieDataPath = rReplace(obj.movieDataPath_,'/|\',filesep);
            
            %Remove ending file separators
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            oldMovieDataPath = regexprep(oldMovieDataPath,endingFilesepToken,'');
            newMovieDataPath = regexprep(newMovieDataPath,endingFilesepToken,'');
            
            %Compare old and new movie paths to detect common tree
            maxNumEl=min(numel(oldMovieDataPath),numel(newMovieDataPath));
            strComp = (oldMovieDataPath(end:-1:end-maxNumEl+1)==newMovieDataPath(end:-1:end-maxNumEl+1));

            %Extract the old and new root directories
            sizeCommonBranch=find(~strComp,1); 
            if isempty(sizeCommonBranch), sizeCommonBranch=maxNumEl+1; end
            oldRootDir=obj.movieDataPath_(1:end-sizeCommonBranch+1);
            newRootDir=newMovieDataPath(1:end-sizeCommonBranch+1);

            % Relocate movie data and channel paths
            for i=1:numel(obj.channels_),
                obj.channels_(i).channelPath_=relocatePath(obj.channels_(i).channelPath_,oldRootDir,newRootDir);
            end
            obj.outputDirectory_=relocatePath(obj.outputDirectory_,oldRootDir,newRootDir);
            obj.movieDataPath_=newMovieDataPath;
            
            % Relocate paths in processes input/output as well as function 
            % and visual parameters 
            for nProc=1:numel(obj.processes_), 
                newInFilePaths = relocatePath(obj.processes_{nProc}.inFilePaths_,oldRootDir,newRootDir);
                obj.processes_{nProc}.setInFilePaths(newInFilePaths);
                newOutFilePaths = relocatePath(obj.processes_{nProc}.outFilePaths_,oldRootDir,newRootDir);
                obj.processes_{nProc}.setOutFilePaths(newOutFilePaths);
                newFunParams = relocatePath(obj.processes_{nProc}.funParams_,oldRootDir,newRootDir);
                obj.processes_{nProc}.setPara(newFunParams);
                newVisualParams = relocatePath(obj.processes_{nProc}.visualParams_,oldRootDir,newRootDir);
                obj.processes_{nProc}.setVisualParams(newVisualParams);
            end
            
        end
        

        function fileNames = getImageFileNames(obj,iChan)
            % Retrieve the names of the images in a specific channel
            if nargin < 2 || isempty(iChan)
                iChan = 1:numel(obj.channels_);
            end            
            if isnumeric(iChan) && min(iChan)>0 && max(iChan) <= ...
                    numel(obj.channels_) && isequal(round(iChan),iChan)                
                fileNames = arrayfun(@(x)(imDir(obj.channels_(iChan(x)).channelPath_)),1:numel(iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.nFrames_)                    
                    error('Incorrect number of images found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
            
        end
        
        function chanPaths = getChannelPaths(obj,iChan)
            %Returns the directories for the selected channels as a cell
            %array or character strings
            nChanTot = numel(obj.channels_);
            if nargin < 2 || isempty(iChan)                
                iChan = 1:nChanTot;
            elseif (max(iChan) > nChanTot) || min(iChan) < 1 || ~isequal(round(iChan),iChan)
                error('Invalid channel index specified! Cannot return path!');                
            end
            nChan = numel(iChan);
            chanPaths = cell(1,nChan);
            for j = 1:nChan                
                %Get the path
                chanPaths{j} = obj.channels_(iChan(j)).channelPath_;                
            end
            
        end
        
        function flag = save(MD)
           % Save the movie data to disk.
           %
           % This function check for the values of the path and filename.
           % If empty, it launches a dialog box asking where to save the
           % MovieData object. If a MAT file already exist, copies this MAT
           % file before saving the movieData. The MovieData variable is
           % saved as 'MD'
           %
           % OUTPUT:
           %    flag - a flag returning the status of the save method
           %
           %
           % Sebastien Besson, 4/2011
           
           % If movieDataPath_ or movieDataFileName_ are empty fields,
           % start a dialog box asking where to save the MovieData
           if isempty(MD.movieDataPath_) || isempty(MD.movieDataFileName_)
               if ~isempty(MD.movieDataPath_),
                   defaultDir=MD.movieDataPath_;
               elseif ~isempty(MD.outputDirectory_)
                   defaultDir=MD.outputDirectory_;
               else
                   defaultDir =pwd;
               end
               [file,path] = uiputfile('*.mat','Find a place to save your movie data',...
                   [defaultDir filesep 'movieData.mat']);
               
               if ~any([file,path]), flag=0; return; end

               % After checking file directory, set directory to movie data
               MD.movieDataPath_=path;
               MD.movieDataFileName_=file;
           end
           
           %First, check if there is an existing file. If so, save a
           %backup. Then save the MovieData as MD
           movieDataFullPath = [MD.movieDataPath_ filesep MD.movieDataFileName_];     
           if exist(movieDataFullPath,'file');
               copyfile(movieDataFullPath,[movieDataFullPath(1:end-3) 'old']);
           end

           save(movieDataFullPath,'MD');
           flag=1;
        end
        
        function setFig = edit(obj)
            setFig = movieDataGUI(obj);
        end

    end
    methods(Static)
        function checkValue=checkValue(property,value)
            % Test the validity of a property value
            %
            % Declared as a static method so that the method can be called
            % without instantiating a MovieData object. Should be called
            % by any generic set method.
            %
            % INPUT:
            %    property - a Channel property name (string)
            %
            %    value - the property value to be checked
            %
            % OUTPUT:
            %    checkValue - a boolean containing the result of the test
            %
            % Sebastien Besson, 4/2011
            
            if iscell(property) 
                checkValue=cellfun(@(x,y) MovieData.checkValue(x,y),property,value);
                return
            end
            
            switch property
                case {'channels_'}
                    checkTest=@(x) isa(x,'Channel');
                case {'outputDirectory_','movieDataPath_','movieDataFileName_','notes_'}
                    checkTest=@(x) ischar(x);
                case {'pixelSize_', 'timeInterval_','numAperture_','magnification_','binning_'}
                    checkTest=@(x) all(isnumeric(x)) && all(x>0);
                case {'camBitdepth_'}
                    checkTest=@(x) isnumeric(x) && ~mod(x, 2);
                case {'owner_'}
                    checkTest= @(x) isa(x,'MovieData');
            end
            checkValue = isempty(value) || checkTest(value);
        end
    end
                
end