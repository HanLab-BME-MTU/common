function status = checkMovieProcedure(movieData,procedure,channels)

% 
% status = checkMovieStatus(movieData,procedure,channels)
% 
% This function checks the status of the requested procedure. If a
% specialized function exists to check this, it is called. If not only the
% basics are checked: Status flag, directory existance, correct field
% format, and file existence if field fileName is specified.
%
% If the specified procedure was completed successfully, status is returned
% as true, and false otherwise. If the named procedure has not been
% correctly specified (as an identically named field) in the movieData,
% false is returned.
% 
% If a function exists named "checkMovie[procedure].m" it will be called
% and its result returned. In this case, this function will do no checking
% of its own.
%
% 
% Input:
% 
%   movieData - The structure describing the movie, as created with
%   setupMovieData.m
%
%   procedure - A character string. The name of the procedure. This should
%   also be the name of a field in the movieData, and a directory within
%   the movie's analysis directory.
%
%   channels - The specific channels of the movie to check the procedure
%   for. Only applies to some procedures. (Optional)
% 
% Hunter Elliott
% 11/2009
% (I know, I should have written this a long time ago....)  

if nargin < 2 || isempty(movieData) || isempty(movieData)
    error('Must input a movieData and specify a procedure to check!')
end
if ~ischar(procedure)
    error('The procedure to check must be specified as a character string!')
end
if nargin < 3
    channels = [];
end

status = false;

%Check if a specialized function exists to check this 
fName = ['checkMovie' upper(procedure(1)) procedure(2:end)]; %The function name should have the procedure capitalized
if exist(fName) == 2        %#ok<EXIST> No second input argument for checking function existence!!!!
    %use this function to check the status
    if nargin(fName) == 2 && ~isempty(channels)    
        
        %If specific channels were requested
        eval(['status = ' fName '(movieData,channels);'])
        return
        
    else
        
        eval(['status = ' fName '(movieData);'])
        return
        
    end
else
    %If no special function, just check the basics. 

    procedure = [lower(procedure(1)) procedure(2:end)]; %In case they capitalized it    
    
    if isfield(movieData,procedure) && isfield(movieData.(procedure),'status') ...
        && movieData.(procedure).status == 1 ...
        && isfield(movieData.(procedure),'directory') && ...
        exist([movieData.(procedure).directory],'dir')
    
        %If a file is specified, make sure it exists
        if isfield(movieData.(procedure),'fileName') && ...
                ~exist([movieData.(procedure).directory filesep movieData.(procedure).fileName],'file')
            status = false;
        else
            status = true;
        end
    end    
end
    


