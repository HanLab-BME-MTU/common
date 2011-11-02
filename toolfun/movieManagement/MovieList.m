classdef MovieList < MovieObject
    % MovieList
    % A class to handle a list of MovieData objects
    
    properties (SetAccess = protected, GetAccess = public)
        
        movieDataFile_       % Cell array of movie data's directory
        movieListPath_       % The path where the movie list is saved
        movieListFileName_   % The name under which the movie list is saved
        
    end
    properties(Transient = true);
        movies_              % Cell array of movies
    end
    
    methods
        function obj = MovieList(movies,outputDirectory, varargin)
            % Constructor for the MovieList object
            
            if nargin > 0
                if iscellstr(movies)
                    if size(movies, 2) >= size(movies, 1)
                        obj.movieDataFile_ = movies;
                    else
                        obj.movieDataFile_ = movies';
                    end
                elseif isa(movies, 'MovieData')
                    obj.movieDataFile_ = arrayfun(@(x) [x.getPath filesep x.getFilename],...
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
        
        
        %% Set/get methods
        function path = getPath(obj)
            path = obj.movieListPath_;
        end
        
        function setPath(obj, value)
            obj.movieListPath_ = value;
        end
        
        function set.movieListPath_(obj, path)
            % Format the path
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            path = regexprep(path,endingFilesepToken,'');
            obj.checkPropertyValue('movieListPath_',path);
            obj.movieListPath_ = path;
        end
        
        function path = getFilename(obj)
            path = obj.movieListFileName_;
        end
        
        function setFilename(obj, filename)
            obj.movieListFileName_ = filename;
        end
        
        function set.movieListFileName_(obj, filename)
            obj.checkPropertyValue('movieListFileName_',filename);
            obj.movieListFileName_ = filename;
        end
        
        
        
        %% Set/get methods
        function movies = getMovies(obj,index)
            % Retrieve the movies from a movie list
            
            if nargin<2 || isempty(index), index = 1:numel(obj.movieDataFile_); end
            movies = cell(numel(index),1);
            for i=index
                movies{i} = MovieData.load(obj.movieDataFile_{i});
            end
        end
        
        
        
        %% Sanitycheck/relocation
        function movieException = sanityCheck(obj, path, filename,askUser)
            % Sanity Check: (Exception 1 - 4)   throws EXCEPTION!
            %
            % ML.sanityCheck
            % ML.sanityCheck(movieListPath, movieListFileName)
            %
            % Assignments:
            %       movieListPath_
            %       movieListFileName_
            %
            % Output:
            %       movieException - cell array of exceptions corresponding to
            %       user index
            
            
            % Ask user by default for relocation
            if nargin < 4, askUser = true; end
            % Calls the superclass sanityCheck (for relocation)
            if nargin>1
                askUser = sanityCheck@MovieObject(obj, path, filename,askUser);
            end
            
            % Apply sanityCheck to all components
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
            
            % Concatenate and throw exceptions if movie loading failed
            if ~all(cellfun(@isempty,movieException)),
                ME = MException('lccb:ml:sanitycheck','Failed to load movie(s)');
                for i=find(~cellfun(@isempty,movieException));
                    ME = ME.addCause(movieException{i});
                end
                throw(ME);
            end
            
            % Save object
            obj.save();
        end
        
        function relocate(obj,newPath)
            % Call superclass relocate function
            [oldRootDir newRootDir]=relocate@MovieObject(obj,newPath);
            
            % Update movie paths
            for i=1:numel(obj.movieDataFile_);
                obj.movieDataFile_{i} = relocatePath(obj.movieDataFile_{i},oldRootDir,newRootDir);
            end
        end
        
        %SB: I am not sure we want to modify the movies list dynamically
        %now we associate packages/processes to a given list (cf MovieData)
%         function removeMovieDataFile(obj, index)
%             % Input:
%             %    index - the index of moviedata to remove from list
%             l = length(obj.movieDataFile_);
%             if any(arrayfun(@(x)(x>l), index, 'UniformOutput', true))
%                 error('User-defined: Index exceeds the length of movie data file.')
%             else
%                 obj.movieDataFile_(index) = [];
%             end
%         end
%         
%         function addMovieDataFile (obj, movie)
%             % Input:
%             %    movie - an array of MovieData objects
%             %            an array of MovieList objects
%             
%             assert( ~iscell(movie), 'User-defined: input cannot be a cell array. It should be a MovieData or MovieList object array.')
%             
%             % Check input data type
%             temp = arrayfun(@(x)(isa(x, 'MovieData')||isa(x, 'MovieList')), movie, 'UniformOutput', true);
%             assert( all(temp), 'User-defined: Input should be a MovieData or MovieList object array')
%             
%             % If no duplicate, add movie data path to MovieList object
%             if isa(movie(1), 'MovieData')
%                 
%                 for i = 1:length(movie)
%                     
%                     if ~any(strcmp(obj.movieDataFile_, [movie(i).getPath()  movie(i).getFilename()]))
%                         obj.movieDataFile_{end+1} = [movie(i).getPath()  movie(i).getFilename()];
%                     end
%                 end
%                 
%             else % MovieList array
%                 for i = 1:length(movie)
%                     
%                     exist = obj.movieDataFile_;
%                     new = movie(i).movieDataFile_;
%                     % temp(0-1 array): 1 - duplicate, 0 - not duplicate
%                     temp = cellfun(@(z)any(z), cellfun(@(x)strcmp(x, exist), new, 'UniformOutput', false), 'UniformOutput', true);
%                     
%                     obj.movieDataFile_ = [obj.movieDataFile_ movie(i).movieDataFile_];
%                 end
%             end
%         end
        
        
    end
    
    methods(Static)
        
        function status = checkProperty(property)
            % Returns true/false if the non-empty property is writable
            status = checkProperty@MovieObject(property);
            if any(strcmp(property,{'movieListPath_','movieListFileName_'}))
                stack = dbstack;
                if any(cellfun(@(x)strcmp(x,'MovieList.sanityCheck'),{stack.name})),
                    status  = true;
                end
            end
        end
        
        % SB: note outputDirectory_ and notes_ should be abstracted to
        % the MovieObject interface but I haven't found a cleanb way to
        % achieve that.
        function status=checkValue(property,value)
            % Return true/false if the value for a given property is valid
            
            if iscell(property)
                status=cellfun(@(x,y) MovieList.checkValue(x,y),property,value);
                return
            end
            
            switch property
                case {'movieListPath_','movieListFileName_','outputDirectory_','notes_'}
                    checkTest=@ischar;
            end
            status = isempty(value) || checkTest(value);
        end
    end
end