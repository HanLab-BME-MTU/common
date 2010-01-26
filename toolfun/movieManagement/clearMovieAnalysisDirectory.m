function movieData = clearMovieAnalysisDirectory(movieData,processNames,clearMiscFiles)

% 
% movieData = clearMovieAnalysisDirectory(movieData)
% 
% WARNING: BE CAREFUL WITH THIS FUNCTION!!!!
% This function "cleans up" the analysis directory of the input movie(s). It
% deletes any processing/analysis results from the selected process(es).
% 
% 
% Input:
% 
%   movieData - The structure describing the movie to cleanup, as created
%   using setupMovieData.m, OR a cell-array of movieData structures if
%   multiple movies are to be cleared.
% 
%   processNames - A character array with the name of the process to clear,
%   or a cell-array with the names of multiple processes to clear. The
%   process specified should be a field of the movieData structure, and
%   should have a sub-field called "directory" which specifies the location
%   of it's results.
%
%   clearMiscFiles - True/False. If set to true, all files in the movie's
%   analysis directory (other than the movieData.mat file) will be removed.
%   Sub-directories will not be removed, unless specified in processNames
%   (See above).
%   Default is false.
%
%
%
% Output:
%
%   movieData - The cleared movieData, with the fields specified by
%   processNames removed. The actual folders and files from these analysis
%   will also be removed from the movie's analysis directory.
%
%
% Hunter Elliott
% 1/2010
%

%% ------ Input ------ %%

if nargin < 2 || isempty(processNames)
    error('Must specify a movieData and a process name to clear! Please specify the processNames input!')
end

%If single movieData or processes were input, convert to cell array.
if ~iscell(movieData)
    movieData = {movieData};
end

if ~iscell(processNames)
    processNames = {processNames};
end

if nargin < 3 || isempty(clearMiscFiles)
    clearMiscFiles = false;
end


%% ------ Clear all movies/processes ------ %%

nMov = length(movieData);
nProc = length(processNames);


for iMov = 1:nMov
    
    
    for iProc = 1:nProc                 
    
        %Look for the current process in the moviedata
        if isfield(movieData{iMov},processNames{iProc})
            
            %Check for the directory specifier
            if isfield(movieData{iMov}.(processNames{iProc}),'directory') && ...
                    exist(movieData{iMov}.(processNames{iProc}).directory,'dir')
                
                disp(['Clearing process "' processNames{iProc} '" for movie ' num2str(iMov)])
                %Remove the directory
                rmdir(movieData{iMov}.(processNames{iProc}).directory,'s')
                
              
            else
                disp(['Could not find directory for process "' processNames{iProc} '" in movieData ' num2str(iMov) ' - can''t clear directory!'])            
            end
            
            %Remove the field
            movieData{iMov} = rmfield(movieData{iMov},processNames{iProc});
            
            
        else
            disp(['Could not find process "' processNames{iProc} '" in movieData ' num2str(iMov) ' - clearing nothing!'])            
        end
        
        if clearMiscFiles
            
            disp(['Clearing all extraneous files from analysis directory for movie ' num2str(iMov) '!'])
            
            %If requested, clear other miscellaneous files from
            %analysis directory.
                    
            %Find all the files
            allFiles = dir([movieData{iMov}.analysisDirectory filesep '*.*']);%Works in both windows and linux
            allFiles = allFiles(arrayfun(@(x)(~x.isdir && ... %But will return the . and .. directories in linux, so remove them...
                ~(strcmp(x.name,'movieData.mat') || strcmp(x.name,'movieData.old') )),allFiles)); %And remove the movieData.mat and movieData.old from the list, so they aren't deleted
                    
            %Remove them all!
            arrayfun(@(x)(delete([movieData{iMov}.analysisDirectory filesep x.name])),allFiles)
                
        end    
    
        
    end
   
   updateMovieData(movieData{iMov});
    
end


