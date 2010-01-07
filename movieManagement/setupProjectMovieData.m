function setupProjectMovieData(projectFolder,pixSize,timeInterval)


% setupProjectMovieData(projectFolder,pixSize,timeInterval)
%
%Description:
%A shitty little function that takes movie folders (which must be set up as
%done in setMovieFolders.m) and creates their movieData using the same
%pixel size and time interval on each.
%
%ALPHA VERSION use with care
%
%Hunter Elliott 2/2009
%
%Input:
%projectFolder - this is the folder which contains all the movie folders.
%This is the same folder you specified when using setupMovieFolders.m
%pixSize - image pixel size in nanometers for all movies in folder
%timeInterval - time interval between images in seconds for all movies in
%folder



%Get the folders for each movie
movieFolders = dir([projectFolder filesep]);

nFiles = length(movieFolders);
isDir = false(nFiles,1);
for j = 1:nFiles
    if movieFolders(j).isdir && ~strcmp(movieFolders(j).name,'.') && ~strcmp(movieFolders(j).name,'..')
        isDir(j) = true;
    end
end
   

iDir = find(isDir);
nMovies = length(iDir);


%Set common properties
movieData.pixelSize_nm = pixSize;
movieData.timeInterval_s = timeInterval;

%Go through each and try to set up the movieData
for j = iDir'
    
    goAhead = 0;
    
    disp(['Setting up movie ' num2str(j-min(iDir)+1) ' of ' num2str(nMovies)])
    
    %Reset the movie data
    movieData.imageDirectory = [];
    movieData.analysisDirectory = [];
    if isfield(movieData,'channelDirectory');
        movieData = rmfield(movieData,'channelDirectory');
    end        

    %Switch to folder
    cd([projectFolder filesep movieFolders(j).name]);
    
    %Check for an existing movieData
    if ~exist('movieData.mat','file')
    
        %Set this as the analysis folder
        movieData.analysisDirectory = pwd;    

        %Look for images folder, make sure it contains images/channels
        imCheck = dir('*images*');    
        nDir = length(imCheck);
        dirGood = zeros(1,nDir);    

        %Check the folders that were found
        for k = 1:length(imCheck)        
            %If its a folder which contains images
            if imCheck(k).isdir && ~isempty(dir([imCheck(k).name filesep '*.tif']))
                dirGood(k) = 1;
            end                           
        end

        %Now finish setting up the movieData             
        if sum(dirGood) == 1
            movieData.imageDirectory = [movieData.analysisDirectory filesep imCheck(dirGood).name];        
            goAhead = 1;
        elseif sum(dirGood) > 1
            disp(['More than one sub-folder named "images" in ' movieData.analysisDirectory '...  cannot setup movieData!!'])
        elseif sum(dirGood) == 0        

            %Check if it is a multi-channel movie
            isChanDir = zeros(1,nDir);
            for k = 1:nDir                        
                %Check for sub-directories containing images                
                tmpDir = dir([movieData.analysisDirectory filesep imCheck(k).name filesep]);
                nImages = zeros(length(tmpDir),1);
                for m = 1:length(tmpDir)
                    if tmpDir(m).isdir && ~strcmp(tmpDir(m).name,'.') && ~strcmp(tmpDir(m).name,'..')
                        nImages(m,k) = length(dir([movieData.analysisDirectory filesep imCheck(k).name filesep tmpDir(m).name filesep '*.tif']));
                    end
                end

                %If all the sub-directories contain same number of images
                if any(nImages(:,k)>0)
                    isChanDir(k) = 1;
                end                                    
            end

            %Make sure there was only one found
            if sum(isChanDir) == 1

                %Set up the folder and sub-folders
                movieData.imageDirectory = [movieData.analysisDirectory filesep imCheck(isChanDir).name];
                tmpDir = dir([movieData.imageDirectory filesep]);
                nGood = 0;
                for m = 1:length(tmpDir)
                    if nImages(m,isChanDir) > 0
                        nGood = nGood + 1;
                        movieData.channelDirectory{nGood} = tmpDir(m).name;
                    end
                end
                goAhead = 1;

            elseif sum(isChanDir) > 1
                disp(['More than one images folder with sub-directories containing images was found in ' movieData.analysisDirectory '...  cannot setup movieData!!'])
            else
                disp(['Could not find a suitable image directory in ' movieData.analysisDirectory '...  cannot setup movieData!!'])                        
            end




            %Check if a frame number is specified in the directory name. This
            %is basically only for use with stimulation/drug perfusion
            %experiments. 



            %Get only the directory from the path:
            r1 = regexpi(movieData.analysisDirectory,filesep);                
            iLastFileSep = max(r1);
            fName = movieData.analysisDirectory(iLastFileSep+1:end);
            r2 = regexpi(fName,'_');
            iLastUscore = max(r2);        
            termString = fName(iLastUscore+1:end);

            %If there was something after the underscore...
            if ~isempty(termString)

                frameNum = str2num(termString);                        

                %Check if it is an integer, and if so store it in movieData
                if ~isempty(frameNum) && ~isnan(frameNum) && round(frameNum) == frameNum
                    movieData.stimulation.iFrame = frameNum;
                end        

            end


        end

        %Specify the shade correction directories
        movieData.FRETprocessing.fretShadeDirectory = [movieData.analysisDirectory filesep 'fretShade'];
        movieData.FRETprocessing.dicShadeDirectory = [movieData.analysisDirectory filesep 'dicShade'];    
        movieData.FRETprocessing.donorShadeDirectory = [movieData.analysisDirectory filesep 'cfpShade'];

        if goAhead
            %Finalize and save the movieData
            setupMovieData(movieData);
        end    

    else
        disp(['Movie ' num2str(j) ' of ' num2str(nMovies) ' already had a movieData.mat! Doing nothing...'])
        
    end
end
%Return to original folder
cd(projectFolder);