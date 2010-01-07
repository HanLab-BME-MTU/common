function movieData = cropMovie(movieData,cropPoly,cropChannels,viewChannels,parentDir)

%
% movieData = cropMovie(movieData, cropPoly)
% 
% This creates a ROI in theinput movie using the input polygon to create a
% mask. If no polygon is input, the user is asked to create one by clicking
% with the mouse
% 
% 
% Input: 
%
%   cropPoly - The polygon to use to crop  the movie. Optional. If not
%   input, the user is asked to create. 
% 
% 
%   cropChannels - The integer index of the channels to crop.
% 
%   viewChannels - The integer index of the channel to view for interactive
%   cropping
%   
%   parentDir   -   The directory to write the cropped images and the new
%   movieData to.
% 
% %Hunter Elliott, 3/2009
% 





%% ------ Parameters ------ %%%

nFrames = 3;%Number of frames to preview for cropping


%% ----- Input ----- %%

if nargin < 1
    %If no movieData is input, the user is asked to select an analysis
    %directory and an image directory below when the moviedata is setup
   movieData = [];
end

if nargin < 3 || isempty(cropChannels)
    %If no channels were input, ask the user which ones to crop    
    [cropChannels,OK] = listdlg('PromptString','Which channel(s) do you want to crop?',...
            'SelectionMode','multiple','ListString',movieData.channelDirectory);
    if ~OK
        return
    end
end

if nargin < 4 || isempty(viewChannels)
    viewChannels = cropChannels(1); %Default is to only view first channel
end

if nargin < 5
    parentDir = [];
end

%Set up and verify the movie data.
movieData = setupMovieData(movieData,'ROI');

nCrop = length(cropChannels);
nView = length(viewChannels);

%If no polygon was input, ask the user to click one
if nargin < 2 || isempty(cropPoly)
    
    retryCrop = true;
    firstTime = true;
    while retryCrop
    
        if firstTime
            %First display images from the movie
            iFrames = round(linspace(1,movieData.nImages(cropChannels(1)),nFrames));
            for j =  1:nFrames                        
                for k = 1:nView
                    subplot(nView,nFrames, (k-1)*nFrames + j )
                    hold on
                    if k == 1 && j == 1
                        title(movieData.analysisDirectory)
                    end

                    axHandle = gca;
                    imageViewer(movieData,'AxesHandle',axHandle,'Frame',iFrames(j)...
                        ,'Channel',movieData.channelDirectory{viewChannels(k)})
                    caxis(caxis ./ 1.8);

                end        
            end
            colormap gray
            firstTime = false;
        end
        %Let the user click a crop outline in the first image of first channel
        subplot(nView,nFrames,1)        
        pHan = impoly(gca);        
        cropPoly = getPosition(pHan);
        delete(pHan);
        %Now draw this polygon on all the frames
        for j =  1:nFrames                        
            for k = 1:nView
                subplot(nView,nFrames, (k-1)*nFrames + j )                        
                fill(cropPoly(:,1),cropPoly(:,2),'r','FaceAlpha',.15,'EdgeColor','r')                                
            end        
        end
    
        %Ask the user if they like their crop or not        
        button = questdlg('Do you like your crop selection?','Confirmation','Yes, Crop Away!','No, Let me try Again.','Abort!','Abort!');
        
        switch button
            
            case 'Yes, Crop Away!'
                retryCrop = false;                            
                
            case 'Abort!'
                return
        end
        
        
    end        
    
end

%% ------ Cropping ----- %%
%Go through each requested channel and each image and crop them

%TEMP - WHAT IF DIR NEEDS TO BE CREATED/NAMED get name from old dirs??

if isempty(parentDir) % if no directory input...
    %Ask the user where they want to put the cropped images
    parentDir = uigetdir(pwd,'Select parent directory for new ROI:');
end
%Make directory for images if necessary
if ~exist([parentDir filesep 'images'],'dir')
    mkdir([parentDir filesep 'images'])
else
    %Delete the old images
    try %TEMP TEMP yeah I know but I'm tired... HLE
        rmdir([parentDir filesep 'images'],'s')
    catch        
    end
    mkdir([parentDir filesep 'images'])
    
end

roiMovieData.imageDirectory = [parentDir filesep 'images'];

wtBar = waitbar(0,'Please wait, cropping images....');

%Go through the images
for iChan = 1:nCrop
    
    %Get current image directory and image names
    imDir = [movieData.imageDirectory filesep movieData.channelDirectory{cropChannels(iChan)} ];
    imNames = dir([imDir filesep '*.tif']);
    nImages = length(imNames); %Allow variable image # per channel
    
    %Store the channel names in the ROI's movieData
    roiMovieData.channelDirectory{cropChannels(iChan)} = movieData.channelDirectory{cropChannels(iChan)};
    
    %Make the channel directory
    mkdir([roiMovieData.imageDirectory filesep roiMovieData.channelDirectory{cropChannels(iChan)}])
    
    for iFrame = 1:nImages
        
        %Load the current image
        currIm = imread([imDir filesep imNames(iFrame).name]);
        
        if iFrame == 1
            [imageM,imageN] = size(currIm);                        
            %Get mask from the polygon
            mask = poly2mask(cropPoly(:,1),cropPoly(:,2),imageM,imageN);
        end
        
        %Crop the image
        currIm(~mask) = 0;
        
        %Remove "extraneous" areas that are outside the crop
        currIm = currIm(max(1,floor(min(cropPoly(:,2) ) ))  : min(ceil(max(cropPoly(:,2) ) ),imageM), ...
                                  max(1,floor(min(cropPoly(:,1) ) ))  : min(ceil(max(cropPoly(:,1) ) ),imageN));
        
        %Write it to the new image directory
        imwrite(currIm,[roiMovieData.imageDirectory filesep roiMovieData.channelDirectory{cropChannels(iChan)} filesep 'crop_' imNames(iFrame).name],'tif')                
        
        
        waitbar( (nImages*(iChan-1) + iFrame ) / (nImages*nCrop),wtBar)
        
    end
    
    
end

close(wtBar);

%% ------- ROI Movie Data Setup ------ %%
%Sets up the movieData for the newly created ROI

%Check if there is a movieData present in the ROI directory, and if not set
%it up.

if ~exist([parentDir filesep 'movieData.mat'],'file')

    %Transfer basic movie info from parent images
    roiMovieData.pixelSize_nm = movieData.pixelSize_nm;
    roiMovieData.timeInterval_s = movieData.timeInterval_s;
    roiMovieData.nImages = movieData.nImages;
    roiMovieData.imageDirectory = [parentDir filesep 'images'];

    if isfield(movieData,'stimulation')
        roiMovieData.stimulation = movieData.stimulation;
    end

    %TEMP - assume that parent is analysis dir
    roiMovieData.analysisDirectory = parentDir;

else %Use the old movieData
   tmp = load([parentDir filesep 'movieData.mat']);
   roiMovieData = tmp.movieData;
    
end
%Save the ROI's movieData
roiMovieData = setupMovieData(roiMovieData);

%Record that an ROI was cropped from the movie in the original movie's
%movieData
nExisting = 0;
if isfield(movieData,'ROI') && isfield(movieData.ROI(1),'analysisDirectory')%If it's been cropped before
    nExisting = length(movieData.ROI);
    isOld = false(1,nExisting);
    for j = 1:nExisting
        if isfield(movieData.ROI(j),'analysisDirectory')
            %Check if this is a re-crop of an existing ROI
            isOld(j) = strcmp(roiMovieData.analysisDirectory,movieData.ROI(j).analysisDirectory);                    
        end
    end
end
     
if exist('isOld','var') 
    if sum(isOld) == 1 %if it found one match        
        iRoi = find(isOld);
    elseif sum(isOld) == 0  %If there was no match but there are previous crops
        iRoi = nExisting  + 1;
    elseif sum(isOld) > 1
        error('Duplicate crops specified in movieData - please check!!')
    end
else %If this is the first crop
    iRoi = 1;
end
    

movieData.ROI(iRoi).analysisDirectory = roiMovieData.analysisDirectory;
movieData.ROI(iRoi).imageDirectory = roiMovieData.imageDirectory;
fName = ['crop data ROI ' num2str(iRoi) '.mat'];
movieData.ROI(iRoi).fileName = fName;
%Save the crop info to the ROI directory
save([movieData.ROI(1).directory filesep fName ],'movieData','cropPoly','cropChannels','imageM','imageN')

movieData.ROI(iRoi).status = 1;
movieData.ROI(iRoi).dateTime = datestr(now);

%Save the input movieData
updateMovieData(movieData)





