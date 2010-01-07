function makeMovieArrayMovies(movieArray)

%Goes through the input movie Array and makes certain movies  of each movie
%whatever figure it out

%Hunter Elliott, 3/2009

nMovies = length(movieArray);

%mFig = figure('Position',[200 240 940 620 ]); %Work desktop large
%mFig = figure('Position',[ 1         125        1024 1082]); %Work desktop fullscreen
mFig = figure('Position',[1 35 1024 660 ]); % Laptop fullscreen

for iMov = 1:nMovies
    
    movieArray{iMov} = refreshMovieData(movieArray{iMov}); %%TEMP?    
    
%     disp(['Making mask movie for movieData ' num2str(iMov) ' of ' num2str(nMovies)])
%     
%      try        
%         makeMaskMovie(movieArray{iMov},[1 2],mFig)
%     catch
%         disp(['Error in mask  movie ' num2str(iMov)])        
%     end
    
    disp(['Making activity overlay movie for movieData ' num2str(iMov) ' of ' num2str(nMovies)])
    clf
    
    try        
        makeActivityOverlayMovie(movieArray{iMov},4,[],mFig);        
    catch errMess
        disp(['Error in activity movie ' num2str(iMov) ' ' errMess.message])        
    end       
%    
%     disp(['Making activity movie for movieData ' num2str(iMov) ' of ' num2str(nMovies)])
%     clf
%     
%     try        
%         makeActivityMovie(movieArray{iMov},3,[],mFig);        
%     catch errMess
%         disp(['Error in activity movie ' num2str(iMov) ' ' errMess.message])        
%     end       
%     
    
    
%     figure(mFig)
%     clf
    
    
%     disp(['Making window movie for movieData ' num2str(iMov) ' of ' num2str(nMovies)])
%     try
%         makeWindowTestingMovie(movieArray{iMov},[],[],2)        
%     catch errMess
%         disp(['Error in window movie ' num2str(iMov) ' : ' errMess.message])        
%     end
    
end

close(mFig);

