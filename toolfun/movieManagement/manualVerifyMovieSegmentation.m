function movieData = manualVerifyMovieSegmentation(movieData,iChannels)

% movieData = manualValidateMovieMasks(movieData,iChannels)
% 
% This function shows several plots and movies which allow
% the user to manually verify that the segmentation is correct/as desired.
% 
% Input:
% 
%   movieData - The movieData structure describing the movie.
%
%   iChannels - The indices of the mask channels to verify. Positive
%                  integer. Optional - if not input, the user will be
%                  asked.
%
% Output:
%
%   movieData - The updated movieDat with the user's evluation stored in it
%               in the field movieData.classification.verify
%
% Hunter Elliott, 10/2009
%

%% ----------- Init -------- %%


%Check the movieData
movieData = setupMovieData(movieData);

if nargin < 2
    iChannels = [];
end

%Check that the segmentation has been completed successfully:
if ~checkMovieMasks(movieData,iChannels)
    error('You must run segmentation first! Check movieData or re-run segmentMovie.m!')
end



%% ------- User verification of Segmentation via Mask Movie ------ %%


%Check that the mask movie has been made
if ~checkMovieMaskMovie(movieData)
    disp('Mask movie has not been made yet, making movie...')
    movieData = makeMaskMovie(movieData,iChannels);    
    
    %Check that it succeeded
    if ~checkMovieMaskMovie(movieData)
        error('Must make mask movie with makeMaskMovie.m to continue!')
    end    
end 

%Compare the date/time on the segmentation and the mask movie
if datenum(movieData.classification.dateTime) > datenum(movieData.classification.movie.dateTime)
    bPressed = questdlg('The segmentation is newer than the mask movie!',...
        'User Mask Validation','Re-Make Movie','Ignore','Abort','Abort');
    
    switch bPressed
        
        case 'Re-Make Movie'
            
            movieData = makeMaskMovie(movieData,iChannels);
            
        case 'Abort'
            
            return
    end
end



%  Open the mask movie

%Tell the user what's going on:
mb=msgbox('Please view the mask movie, and then close the viewing application when finished.','modal');
uiwait(mb);

if ispc
    fStat = system([movieData.analysisDirectory filesep movieData.classification.movie.fileName]);        
    
elseif isunix
    cd(movieData.analysisDirectory)
    fStat = system(['totem '  movieData.classification.movie.fileName]); %not ideal, but what else?
end

if fStat ~= 0
    error('Problem displaying mask movie!')
end

%Ask the user how the masks looked
bPressed = questdlg('Were the masks correct?','User Mask Validation','Yes','No','No');

if strcmp(bPressed,'Yes')
    movieData.classification.manualValidateMasks = {'Good'};    
else
    movieData.classification.manualValidateMasks = {'Bad'};
end



%Save the movieData
updateMovieData(movieData);

