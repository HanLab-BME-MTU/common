classdef MovieList < MovieObject
    % Concrete implementation of MovieObject for a list of movies
    
    properties
        movieListFileName_   % The name under which the movie list is saved
        movieListPath_       % The path where the movie list is saved
    end
    properties (SetAccess = protected)
        movieDataFile_       % Cell array of movie data's directory
    end
    properties(Transient = true);
        movies_              % Cell array of movies
    end
    
    methods
        function obj = MovieList(movies,outputDirectory, varargin)
            % Constructor for the MovieList object
            
            if nargin > 0
                if iscellstr(movies)
                    obj.movieDataFile_ = movies(:)';
                elseif iscell(movies) && all(cellfun(@(x)isa(x,'MovieData'),movies))
                    obj.movieDataFile_ = cellfun(@(x) fullfile(x.getPath,x.getFilename),...
                        movies,'UniformOutput',false);
                elseif isa(movies, 'MovieData')
                    obj.movieDataFile_ = arrayfun(@(x) fullfile(x.getPath,x.getFilename),...
                        movies,'UniformOutput',false);
                else
                    error('lccb:ml:constructor','Movies should be a cell array or a array of MovieData');
                end
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
        end
        
        
        %%  Set/get methods

        function set.movieListPath_(obj, path)
            % Set movie list path
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            path = regexprep(path,endingFilesepToken,'');
            obj.checkPropertyValue('movieListPath_',path);
            obj.movieListPath_ = path;
        end
        
        function set.movieListFileName_(obj, filename)
            obj.checkPropertyValue('movieListFileName_',filename);
            obj.movieListFileName_ = filename;
        end
               
        function movies = getMovies(obj,varargin)
            % Get the movies from a movie list
            
            ip =inputParser;
            allIndex = 1:numel(obj.movieDataFile_);
            ip.addOptional('index',allIndex,@(x) all(ismember(x,allIndex)));
            ip.parse(varargin{:});
            index= ip.Results.index;
            
            movies = cell(numel(index),1);
            for i=index, movies{i} = MovieData.load(obj.movieDataFile_{i}); end
        end
            
        %% Sanity check/relocation
        function movieException = sanityCheck(obj, varargin)
            % Check the sanity of the MovieData objects
            %
            % First call the superclass sanityCheck. Then load the individual 
            % movies in the list (runs sanityCheck on each movie).
            % Save the movie list to disk if run successfully.
            
            % Call the superclass sanityCheck
            if nargin>1, 
                askUser = sanityCheck@MovieObject(obj, varargin{:});
            else
                askUser = true;
            end
            
            % Load movie components (run sanityCheck on each of them)
            movieIndex = 1:numel(obj.movieDataFile_);
            movieException = cell(1, numel(movieIndex));
            for i = movieIndex
                try
                    obj.movies_{i}=MovieData.load(obj.movieDataFile_{i},askUser);
                catch ME
                    movieException{i} = ME;
                    continue
                end
            end
            
            % Throw exception if at least one movie failed during loading
            if ~all(cellfun(@isempty,movieException)),
                ME = MException('lccb:ml:sanitycheck','Failed to load movie(s)');
                for i=find(~cellfun(@isempty,movieException));
                    ME = ME.addCause(movieException{i});
                end
                throw(ME);
            end
            
            obj.save();
        end
        
        function relocate(obj,varargin)
            % Relocate the MovieList object
            
            % Run superclass relocate function (movie list path and analysis components)
            [oldRootDir newRootDir]=relocate@MovieObject(obj,varargin{:});
            
            % Relocate movie paths
            for i=1:numel(obj.movieDataFile_);
                obj.movieDataFile_{i} = relocatePath(obj.movieDataFile_{i},oldRootDir,newRootDir);
            end
        end
    end
    
    methods(Static)
        
         function status=checkValue(property,value)
           % Return true/false if the value for a given property is valid
            
           % Parse input
           ip = inputParser;
           ip.addRequired('property',@(x) ischar(x) || iscell(x));
           ip.parse(property);
           if iscell(property)
               ip.addRequired('value',@(x) iscell(x)&&isequal(size(x),size(property)));
               ip.parse(property,value);
               status=cellfun(@(x,y) MovieList.checkValue(x,y),property,value);
               return
           end
           
           % Get validator for single property
           validator=MovieList.getPropertyValidator(property);
           propName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
           assert(~isempty(validator),['No validator defined for property ' propName]);
           
           % Return result of validation
           status = isempty(value) || validator(value);
        end
        
        function validator = getPropertyValidator(property) 
            validator = getPropertyValidator@MovieObject(property);
            if ~isempty(validator), return; end
            if ismember(property, {'movieListPath_','movieListFileName_'})
                validator=@ischar;
            end
        end
        
        function propName = getPathProperty()
            propName = 'movieListPath_';
        end
        function propName = getFilenameProperty()
            propName = 'movieListFileName_';
        end
    end
end