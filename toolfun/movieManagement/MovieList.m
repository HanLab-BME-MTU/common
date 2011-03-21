classdef MovieList < handle
   
    properties (SetAccess = protected, GetAccess = public)
       
        movieDataFile_ = {}  % Cell array of movie data's directory
        
        movieListPath_
        movieListFileName_
        
    end
    
    methods (Access = public)
        function obj = MovieList(movieDataFile, movieListPath, movieListFileName)
            % movieDataFile: a cell array of file directory
            %                an array of MovieData object
            if nargin > 0
                if iscellstr(movieDataFile)
                    if size(movieDataFile, 2) >= size(movieDataFile, 1)
                        obj.movieDataFile_ = movieDataFile;
                    else
                        obj.movieDataFile_ = movieDataFile';
                    end
                    
                elseif isa(movieDataFile, 'MovieData')
                    for i = 1: length(movieDataFile)
                       obj.movieDataFile_{end+1} = [movieDataFile(i).movieDataPath_ filesep movieDataFile(i).movieDataFileName_]; 
                    end
                else
                    error('User-defined: Please provide a cell array of file directory or an array of MovieData object.')
                end
                
                if nargin > 1
                   if isempty(movieListPath) || ischar(movieListPath)
                       obj.movieListPath_ = movieListPath;
                   else
                       error('User-defined: Movie list path should be a string');
                   end
                end
                
                if nargin > 2
                   if isempty(movieListFileName) || ischar(movieListFileName)
                       obj.movieListFileName_ = movieListFileName;
                   else
                       error('User-defined: Movie list file name should be a string');
                   end
                end
                
            else
                error('User-defined: Please provide movie data file path to creat a movie list.')
            end
            
        end
        
        function [movieException, MDList] = sanityCheck(obj, userIndex,movieListPath, movieListFileName)
            %
            % Sanity Check: (Exception 1 - 4)   throws EXCEPTION!
            %
            % ML.sanityCheck
            % ML.sanityCheck('all')
            % ML.sanityCheck(userIndex)
            % ML.sanityCheck(userIndex, movieListPath, movieListFileName)
            %
            % Assignments:
            %       movieListPath_
            %       movieListFileName_
            %
            % Output:
            %       movieException - cell array of exceptions corresponding to
            %       user index
            %       
            %       MDList - cell array of Movie Data objects 
            %
            
            if nargin < 2
                userIndex = 'all';
            end
            
            if strcmp(userIndex, 'all')
                index = 1:length(obj.movieDataFile_);

            elseif max(userIndex) <= length(obj.movieDataFile_)
                index = userIndex;
            else
                error('User-defined: user index exceed the length of movie data list.')                
            end
            
            movieException = cell(1, length(index));
            MDList = cell(1, length(index));
            askUser=true;

            if nargin > 2  
            % Check if the path and filename stored in the movieList are the same
            % as the ones provided in argument. They can differ if the movieList
            % MAT file has been renamed, move or copy to another location.

                endingFilesepToken = [regexptranslate('escape',filesep) '$'];
                oldPath = regexprep(obj.movieListPath_,endingFilesepToken,'');
                newPath = regexprep(movieListPath,endingFilesepToken,'');
                if  ~strcmp(oldPath, newPath)
                    relocateMsg=sprintf(['The movie list located in \n%s\n has been relocated to \n%s.\n'...
                        'Should I try to relocate the movie data in the list as well?'],oldPath,newPath);
                    confirmRelocate = questdlg(relocateMsg,'Movie Data','Yes to all','Yes','No','Yes');
                    if strcmp(confirmRelocate(1:3),'Yes')
                        obj.relocateMovieList(movieListPath); 
                        askUser = (strcmp(confirmRelocate,'Yes'));
                    else
                        obj.setMovieListPath(movieListPath);
                    end
                end
            
                if  ~strcmp(obj.movieListFileName_, movieListFileName)
                    obj.movieListFileName_ = movieListFileName; 
                end
            end
            

            for i = 1 : length(index)
                
            % Exception 1: MovieData file does not exist  
            
                if ~exist(obj.movieDataFile_{index(i)}, 'file')
                    movieException{i} = MException('lccb:ml:nofile', 'File does not exist.');
                    continue
                end
                
            % Exception 2: Fail to open .mat file                
                try
                    pre = whos('-file', obj.movieDataFile_{index(i)});
                catch ME
                    movieException{i} = MException('lccb:ml:notopen', 'Fail to open file. Make sure it is a MAT file.');
                    continue
                end
                
            % Exception 3: No MovieData object in .mat file
            
                structMD = pre( logical(strcmp({pre(:).class},'MovieData')) );
                switch length(structMD)
                    case 0
                        movieException{i} = MException('lccb:ml:nomoviedata', ...
                            'No movie data is found in selected MAT file.');
                        continue                        
                        
                    case 1
                        data = load(obj.movieDataFile_{index(i)}, '-mat', structMD.name);
                        MDList{i}=data.(structMD.name);
                        
                        
            % Exception 4: More than one MovieData objects in .mat file
            
                    otherwise
                        movieException{i} = MException('lccb:ml:morethanonemoviedata', ...
                            'More than one movie data are found in selected MAT file.');
                        continue     
                end
                
            % Exception 5: Movie Data Sanity Check 
            
                try 
                    [path filename ext]=fileparts(obj.movieDataFile_{index(i)});
                    MDList{i}.sanityCheck(path,[filename ext],askUser)
                    
                catch ME
                    movieException{i} = MException('lccb:ml:sanitycheck', ME.message);
                    continue                     
                    
                end
            end
            
        end
        
        function removeMovieDataFile (obj, index)
            % Input:
            %    index - the index of moviedata to remove from list
            l = length(obj.movieDataFile_);
            if any(arrayfun(@(x)(x>l), index, 'UniformOutput', true))
                error('User-defined: Index exceeds the length of movie data file.')
            else
                obj.movieDataFile_(index) = [];
            end
        end
        
        function editMovieDataFile(obj, index, movieDataFile)
            % Assign the index(i) th movie data path to movieDataFile{i}
            %
            % Input:
            %    index - array of the index of movie data to edit
            %    movieDataFile - cell array of new movie data path
            
            if iscellstr(movieDataFile)
                if size(movieDataFile, 2) < size(movieDataFile, 1)
                    movieDataFile = movieDataFile';
                end
            else
                error('User-defined: input movieDataFile should be a string cell array')
            end
            
            l = length(obj.movieDataFile_);
            if any( arrayfun(@(x)(x > l), index, 'UniformOutput', true) )
               error('User-defined: input index exceeds the length of movie data.') 
            end
            
            assert( length(index) == length(movieDataFile), 'User-defined: the length of input index and movieDataFile must be equal;')
            
            % Assign movie data path 
            obj.movieDataFile_(index) = movieDataFile;
            
        end
        
        function addMovieDataFile (obj, movie)
            % Input:
            %    movie - an array of MovieData objects
            %            an array of MovieList objects
            
            assert( ~iscell(movie), 'User-defined: input cannot be a cell array. It should be a MovieData or MovieList object array.')
            
            % Check input data type
            temp = arrayfun(@(x)(isa(x, 'MovieData')||isa(x, 'MovieList')), movie, 'UniformOutput', true);
            assert( all(temp), 'User-defined: Input should be a MovieData or MovieList object array')
            
            % If no duplicate, add movie data path to MovieList object
            if isa(movie(1), 'MovieData')
                
                for i = 1:length(movie)
                    
                    if ~any(strcmp(obj.movieDataFile_, [movie(i).movieDataPath_  movie(i).movieDataFileName_]))
                        obj.movieDataFile_{end+1} = [movie(i).movieDataPath_  movie(i).movieDataFileName_];
                    end
                end
                
            else % MovieList array
                for i = 1:length(movie)
                    
                    exist = obj.movieDataFile_;
                    new = movie(i).movieDataFile_;
                    % temp(0-1 array): 1 - duplicate, 0 - not duplicate
                    temp = cellfun(@(z)any(z), cellfun(@(x)strcmp(x, exist), new, 'UniformOutput', false), 'UniformOutput', true);
                    
                    obj.movieDataFile_ = [obj.movieDataFile_ movie(i).movieDataFile_];
                end
            end
        end

        function setMovieListPath(obj, path)
            % Set the path to the MAT file containing the movie list
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            obj.movieListPath_ = regexprep(path,endingFilesepToken,''); 
        end
        
        
        function relocateMovieList(obj,newMovieListPath)
            % Relocate all movie data files in movie list if applicable
            
            %Convert temporarily all paths using the local file separator.
            %Remove ending separators 
            oldMovieListPath = rReplace(obj.movieListPath_,'/|\',filesep);
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            oldMovieListPath = regexprep(oldMovieListPath,endingFilesepToken,'');
            newMovieListPath = regexprep(newMovieListPath,endingFilesepToken,'');
            
            %Compare old and new movie paths, detect common tree and
            %extract the old and new root directories
            maxNumEl=min(numel(oldMovieListPath),numel(newMovieListPath));
            strComp = (oldMovieListPath(end:-1:end-maxNumEl+1)==newMovieListPath(end:-1:end-maxNumEl+1));
            sizeCommonBranch=find(~strComp,1); 
            oldRootDir=obj.movieListPath_(1:end-sizeCommonBranch+1);
            newRootDir=newMovieListPath(1:end-sizeCommonBranch+1);

            %Guess paths of movie data files using the old and new root
            newMovieDataPaths = cellfun(@(x) relocatePath(x,oldRootDir,newRootDir),obj.movieDataFile_,'Unif',false);
            changedMovieDataPaths=(~cellfun(@isempty,newMovieDataPaths));           

%             confirmRelocate = questdlg('The location of some movie data has changed. Should I replace the locations of these elements?',...
%                  'Movie List Relocate','Yes','No','Yes');

%             if strcmp(confirmRelocate,'Yes')
            obj.editMovieDataFile(changedMovieDataPaths,newMovieDataPaths)
%             end

           
            obj.setMovieListPath(newMovieListPath);
        end
            
        function saveMovieList(ML)
           % Save movie list to disk. Save a backup of existing copy if
           % applicable.
           
           if exist([ML.movieListPath_ filesep ML.movieListFileName_],'file');
               copyfile([ML.movieListPath_ filesep ML.movieListFileName_],...
                    [ML.movieListPath_ filesep ML.movieListFileName_(1:end-3) 'old'])
           end
          
           save([ML.movieListPath_ filesep ML.movieListFileName_],'ML')
        end
        
    end
    
end

