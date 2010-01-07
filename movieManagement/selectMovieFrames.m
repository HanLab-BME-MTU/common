function movieData = selectMovieFrames(movieData,visChan)

%{
 movieData = selectMovieFrames(movieData or movieArray,visChans)
Allows the user to select which frames to use from the input movie / movie array


visChans - The indices of the channels the user would like to view for frame selection

%Hunter Elliott, 4/2009
%}


%% ------ Input ----- %%

wasSingle = false;

if nargin < 2 || isempty(visChan)
    visChan = 1;
elseif length(visChan) ~= 1|| round(abs(visChan)) ~= visChan
    errordlg('Requested visualization channel must be an integer!',mfilename)
    return
end

if ~iscell(movieData) %If a single moviedata was input, convert to cell array
    movieData = {movieData};
    wasSingle = true;
end
       


nMovies = length(movieData);



%% ---- Frame Selections ---- %%
%Go through the movie(s) and allow the user to select frames

if nMovies > 1
    progressbar(0);
end

for iMov = 1:nMovies
    
    %Check / convert this movieData
    movieData{iMov} = setupMovieData(movieData{iMov});
    
    movieData{iMov}.selectedFrames.status = 0;
    
    disp(['Movie ' num2str(iMov) ' of ' num2str(nMovies) ' :'])
    disp(movieData{iMov}.analysisDirectory)
    
    %open the window selection movie if one exists
    try
        system([movieData{iMov}.analysisDirectory filesep 'windowTestingMovie.mov']);
    catch errMess
        disp(['Could not open window testing movie for ' movieData{iMov}.analysisDirectory ' : ' errMess.message])
    end
    
    
    %Get image file names
    imDir = [movieData{iMov}.imageDirectory filesep movieData{iMov}.channelDirectory{visChan} ];
    imNames = dir([imDir filesep '*.tif']);
    
    
    
    if length(imNames) ~= movieData{iMov}.nImages
        errordlg('Incorrect number of images in selected image channel! Check movie data and image directories!',mfilename)
    end    
    
    %Open the requested channel for viewing using default image viewer    
    system([imDir filesep imNames(1).name]);
    
    
        
    %Now ask the user which frames to keep
    selFrames = listSelectGUI(1:movieData{iMov}.nImages,[],'move');
    
    movieData{iMov}.selectedFrames.iFrames = selFrames;
    movieData{iMov}.selectedFrames.status = 1;
    movieData{iMov}.selectedFrames.dateTime = datestr(now);
    
    updateMovieData(movieData{iMov})
    
    if nMovies > 1
        progressbar(iMov/nMovies);
    end
    
end


