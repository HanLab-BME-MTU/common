function splitProjectMovieImages(projectDir)

%SPLITMOVEIIMAGES splits multi-image .tif files into single tifs and places them in individual directories 
% 
% splitMovieImages(projectDir)
% 
% This function goes through every sub-directory of the directory projectDir
% and if that directory contains .tif files, a sub-directory called
% "Images" is created. Each .tif is the placed in a seperate sub-directory
% of this Images directory which is named after the .tif, and if it is a
% multi-page .tif file the pages are split into seperate, sequentially
% named files.
% 
% Input:
% 
%   projectDir - Parent directory containing sub-directories, each of
%   whichc contains images.
% 
% 
% 
% 
% Hunter Elliott
% 2/2010



fExt = 'tif';%File extension to look for images
storeName = 'originalStacks'; %The name of the directory to put the original stacks in after they have been split



allSub = dir(projectDir); %Check contents of project directory

allSub = allSub([allSub.isdir]' &  ...
    arrayfun(@(x)(~any(strcmp(x.name,{'.','..'}))),allSub)); %Remove all non-directories, and the . & .. directories. This works on PC & linux.

nSub = length(allSub);

for i = 1:nSub
    
    %Check for images in this directory
    im = dir([projectDir filesep allSub(i).name filesep '*.' fExt ]);
    
    if ~isempty(im)
        
        disp(['Splitting images for folder ' num2str(i)])
        
        %Set up the image directory
        imDir = [projectDir filesep allSub(i).name filesep 'images'];        
        mkdir(imDir)
        
        %Make a folder for storing the old stacks
        mkdir([projectDir filesep allSub(i).name filesep storeName]);
        
        %Write the images to seperate sub-dirs of this image dir
        for j = 1:length(im)
            
            %Make the image directory
            mkdir([imDir filesep im(j).name(1:end-4)])%Name the directory after the stack, removing the file extension
            
            %Load all the images
            currIm = stackRead([projectDir filesep allSub(i).name filesep im(j).name]);
            nIm = size(currIm,3); %Check number of images
            %Get number of digits for writing file names
            nDig = floor(log10(nIm)+1);
            %Make the string for formatting
            fString = strcat('%0',num2str(nDig),'.f');
            %Write them all to new dir
            disp(['Splitting "' im(j).name '" into ' num2str(nIm) ' seperate files...'])
            for k = 1:nIm
                imwrite(squeeze(currIm(:,:,k)),[imDir filesep im(j).name(1:end-4) filesep im(j).name(1:end-4) '_' num2str(k,fString) '.' fExt]);
            end
            
            movefile([projectDir filesep allSub(i).name filesep im(j).name],[projectDir filesep allSub(i).name filesep storeName])
            
        end
    else
        disp(['No images found in sub-folder ' allSub(i).name '!!']);
    end
    
    
    
    
end

