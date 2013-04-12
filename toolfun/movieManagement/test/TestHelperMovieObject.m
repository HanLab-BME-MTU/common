classdef TestHelperMovieObject < handle
    methods (Static)
        %% Set up and tear down methods
        function movieList = setUpMovieList(path,movie)
            if nargin<2, movie=TestHelperMovieObject.setUpMovie(path);end
            if ~exist(path,'dir'), mkdir(path); end
            movieList = MovieList(movie,path);
            movieList.setPath(path);
            movieList.setFilename('movieList.mat');
            movieList.sanityCheck;
        end
        
        function movie = setUpMovie(path,channel)
            if nargin<2, channel=TestHelperMovieObject.setUpChannel(path);end
            if ~exist(path,'dir'), mkdir(path); end
            movie = MovieData(channel,path);
            movie.setPath(path);
            movie.setFilename('movieData.mat');
        end
        
        function channel = setUpChannel(path,imSize,nFrames, format)
            if nargin<4, format = 'double'; end
            if nargin<3, nFrames=1;end
            if nargin<2,imSize=[100 200]; end
            if ~exist(path,'dir'), mkdir(path); end
            for i=1:nFrames
                imwrite(zeros(imSize, format),fullfile(path,['test_' num2str(i) '.tif']));
            end
            channel=Channel(path);
        end
        
        function movie = setUpBFMovie(path, imSize, nChan, nFrames, format)
            if nargin<5, format = 'double'; end
            if ~exist(path,'dir'), mkdir(path); end
            bfsave(zeros(imSize(1), imSize(2), nChan, nFrames, format),...
                fullfile(path, 'test.ome.tiff'), 'XYCTZ');
            movie = bfImport(fullfile(path, 'test.ome.tiff'));
        end
        
        function testGetProcessIndex(movieObject)
            allProc= TestHelperMovieObject.getConcreteSubClasses('Process');
            procClass = allProc{1};
            assertEqual(movieObject.getProcessIndex(procClass,1,false),[]);
            
            procConstr = str2func(procClass);
            movieObject.addProcess(procConstr(movieObject,''));
            assertEqual(movieObject.getProcessIndex(procClass,1,false),1);
            
            movieObject.addProcess(procConstr(movieObject,''));
            assertEqual(movieObject.getProcessIndex(procClass,Inf,false),1:2);
            assertEqual(movieObject.getProcessIndex(procClass,2,false),1:2);
        end
        
        
        function testProcessCreation(movieObject)
            concreteProc= TestHelperMovieObject.getConcreteSubClasses('Process');
            for i=1:numel(concreteProc)
                procConstr = str2func(concreteProc{i});
                switch concreteProc{i}
                    case 'SignalProcessingProcess'
                        assertExceptionThrown(@()procConstr(movieObject),'lccb:getDefaultParams:noWindowing')
                    otherwise
                        % Construct new process
                        newprocess = procConstr(movieObject);
                        assertTrue(isa(newprocess,concreteProc{i}));
                        
                        % Add process
                        movieObject.addProcess(newprocess);
                        assertTrue(isa(movieObject.processes_{1},concreteProc{i}));
                        movieObject.deleteProcess(1);
                end
                movieObject.reset
                assertTrue(isempty(movieObject.processes_));
            end
            
        end
        
        
        function testPackageCreation(movieObject)
            % Test full package creation
            
            concretePackages = TestHelperMovieObject.getConcreteSubClasses('Package');
            for i=1:numel(concretePackages)
                % Test package construction
                packageConstr = str2func(concretePackages{i});
                if strcmp(concretePackages{i},'IntegratorPackage') && ...
                        isa(movieObject,'MovieData')
                    assertExceptionThrown(@()procConstr(movieObject),'MATLAB:InputParser:ArgumentFailedValidation')
                else
                    newpackage = packageConstr(movieObject,'');
                    assertTrue(isa(newpackage,concretePackages{i}));
                    
                    % Test package addition
                    movieObject.addPackage(newpackage);
                    assertTrue(isa(movieObject.packages_{end},concretePackages{i}));
                    
                    % Test default process construction
                    crtPackage = movieObject.packages_{1};
                    for j=1:numel(crtPackage.getDefaultProcessConstructors),
                        crtPackage.createDefaultProcess(j);
                        assertTrue(isa(movieObject.processes_{end},...
                            crtPackage.getProcessClassNames{j}));
                    end
                    
                    % Reset movie
                    movieObject.reset
                    assertTrue(isempty(movieObject.processes_));
                    assertTrue(isempty(movieObject.packages_));
                end
            end
        end
        
        function relocatedMoviePath= relocateMovie(movieObject)
            % Copy movie in new location
            relocatedMoviePath = [movieObject.getPath '_relocated'];
            copyfile(movieObject.getPath,relocatedMoviePath);
            rmdir(movieObject.getPath, 's')
        end

        
        function testMultiSetProperties(object,validProperties,validValues)
            set(object,validProperties,validValues);
            set(object,validProperties,validValues);
            
            for i=1:numel(validProperties)
                f= @() set(object,validProperties{i},NaN);
                assertExceptionThrown(f,'lccb:set:readonly');
            end
        end
        

        function concreteSubClasses=getConcreteSubClasses(superclass)
            % List all files in the matlab path ending by Process
            pathList=regexp(path,pathsep,'split');
            f=cellfun(@(x) dir([x filesep '*' superclass '.m']),pathList,'Unif',false);
            f=f(~cellfun(@isempty,f));
            f=vertcat(f{:});
            allClasses = cellfun(@(x) x(1:end-2),{f.name},'Unif',false);
            allClasses = allClasses(cellfun(@(x) ~isempty(meta.class.fromName(x)),allClasses));
            
            % Get concrete subclasses
            subClasses= allClasses(cellfun(@(x) isSubclass(x,superclass),allClasses));
            concreteSubClasses= allClasses(cellfun(@isConcreteClass,subClasses));
            concreteSubClasses= unique(concreteSubClasses);
        end
    end
end