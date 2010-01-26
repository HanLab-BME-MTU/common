function setupMovieFolders(projectFolder)

%Description:
%A quick and dirty script to set up movie data folders assuming it is a CFP
%single chain biosensor with appropriately named images and shade images
%
%INPUTS: The project folder is the super-folder which contains 1 folder for
%each movie with all the images in it (shade, CFP, FRET...)
%
%ALPHA VERSION!! Use with care!!!
%
%Hunter Elliott 2/09
%



%Get the folders for each movie
movieFolders = dir([projectFolder]);
nMovies = length(movieFolders);


%Check isdir on each entry
isDir = false(nMovies,1);
for j = 1:length(movieFolders)
    
    
    if movieFolders(j).isdir && ~strcmp(movieFolders(j).name,'.') && ~strcmp(movieFolders(j).name,'..')    
        isDir(j) = true;
    end

end

goodDir = find(isDir)'; 


%Go through each folder and set up the image folders
for j = goodDir
    
    
    
    disp(['Setting up movie ' num2str(j) ' of ' num2str(nMovies)])
    
    %Switch to movie directory
    cd([projectFolder filesep movieFolders(j).name]);
    
    %Make shade correction folders
    mkdir('cfpShade'),mkdir('dicShade'),mkdir('fretShade')
    %And image folder
    mkdir('images')
    
    disp('making shade folders...');
    
    %Find and move the shade correction images
    %Cfp
    
    disp('moving shade images..')
    
    cfpShadeFiles = dir('*CFPshad*.tif');
    nCfpShade = length(cfpShadeFiles);
    for m = 1:nCfpShade        
        movefile(cfpShadeFiles(m).name,['.' filesep 'cfpShade' filesep]);                
    end
    %fret
    fretShadeFiles = dir('*FRETshad*.tif');
    nfretShade = length(fretShadeFiles);
    for m = 1:nfretShade        
        movefile(fretShadeFiles(m).name,['.' filesep 'fretShade' filesep]);                
    end
    %DIC - sometimes I am sent these... not used though...
    dicShadeFiles = dir('*DICshad*.tif');
    ndicShade = length(dicShadeFiles);
    for m = 1:ndicShade        
        movefile(dicShadeFiles(m).name,['.' filesep 'dicShade' filesep]);                
    end
    
    %Now make the image folders
    cd('images')
    mkdir('CFP'),mkdir('DIC'),mkdir('FRET')
    
    disp('making image folders..')
    
    %And move the image files into them
    
    %Cfp    
    disp('moving image files...')
    cd('..')
    cfpFiles = dir('*CFP*.tif');
    nCfp = length(cfpFiles);
    for m = 1:nCfp        
        movefile(cfpFiles(m).name,['.' filesep 'images' filesep 'CFP' filesep])                
    end
    
    fretFiles = dir('*FRET*.tif');
    nfret = length(fretFiles);
    for m = 1:nfret        
        movefile(fretFiles(m).name,['.' filesep 'images' filesep 'FRET' filesep])                
    end
    
    dicFiles = dir('*DIC*.tif');
    ndic = length(dicFiles);
    for m = 1:ndic        
        movefile(dicFiles(m).name,['.' filesep 'images' filesep 'DIC' filesep])                
    end
    
    %Check if the number of images matches
    
    if ndic == nfret && ndic == nCfp
        disp(['Set up ' num2str(nCfp) ' images successfully'])        
    else
        disp(['Warning: mismatch in number of images!!! nCFP: ' num2str(nCfp) ' nFRET: ' num2str(nfret) ' nDIC: ' num2str(ndic)])
    end
    
    
    
end
