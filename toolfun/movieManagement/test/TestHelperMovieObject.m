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
            movie.sanityCheck;
        end
        
        function channel = setUpChannel(path,imSize,nFrames)
            if nargin<3, nFrames=1;end
            if nargin<2,imSize=[100 200]; end
            if ~exist(path,'dir'), mkdir(path); end
            for i=1:nFrames
                imwrite(zeros(imSize),fullfile(path,['test_' num2str(i) '.tif']));
            end
            channel=Channel(path);
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
            procConstr = cellfun(@str2func,concreteProc,'Unif',false);
            for i=1:numel(procConstr)
                newprocess = procConstr{i}(movieObject);
                assertTrue(isa(newprocess,concreteProc{i}));
                
                movieObject.addProcess(newprocess);
                assertTrue(isa(movieObject.processes_{end},concreteProc{i}));
            end
        end
        
        
        function testPackageCreation(movieObject)
            concretePackages = TestHelperMovieObject.getConcreteSubClasses('Package');
            packageConstr = cellfun(@str2func,concretePackages,'Unif',false);
            for i=1:numel(packageConstr)
                newpackage = packageConstr{i}(movieObject,'');
                assertTrue(isa(newpackage,concretePackages{i}));
                
                
                movieObject.addPackage(newpackage);
                assertTrue(isa(movieObject.packages_{end},concretePackages{i}));
                
                crtPackage = movieObject.packages_{end};
                for j=1:numel(crtPackage.getDefaultProcessConstructors),
                    procConstr = crtPackage.getDefaultProcessConstructors{j};
                    movieObject.addProcess(procConstr(movieObject,crtPackage.outputDirectory_));
                    crtPackage.setProcess(j,movieObject.processes_{end});
                    
                    assertTrue(isa(movieObject.processes_{end},crtPackage.getProcessClassNames{j}));
                end
            end
        end
        
        
        function relocatedMoviePath= relocateMovie(movieObject)
            % Copy movie in new location
            relocatedMoviePath = [movieObject.getPath '_relocated'];
            copyfile(movieObject.getPath,relocatedMoviePath);
        end
        
        function testSetInvalidProperties(object,validProperties)
            for i=1:numel(validProperties)
                f= @() set(object,validProperties{i},NaN);
                assertExceptionThrown(f,'lccb:set:invalid');
                
                f= @() set(object,validProperties{i},0);
                assertExceptionThrown(f,'lccb:set:invalid');
            end
        end
        
        function testMultiSetProperties(object,validProperties,validValues)
            set(object,validProperties,validValues);
            set(object,validProperties,validValues);
            
            for i=1:numel(validProperties)
                f= @() set(object,validProperties{i},NaN);
                assertExceptionThrown(f,'lccb:set:readonly');
            end
        end
        
        
        function testSetMultipleProperties(object,validProperties,validValues)
            set(object,validProperties,validValues);
            for i=1:numel(validProperties)
                assertEqual(object.(validProperties{i}),validValues{i});
            end
        end
        
        
        function testSetIndividualProperties(object,validProperties,validValues)
            for i=1:numel(validProperties)
                set(object,validProperties{i},validValues{i});
                assertEqual(object.(validProperties{i}),validValues{i});
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
        end
    end
end