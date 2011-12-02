function mismatchRgbMerge()
%This function is designed to make RGB images of two channels where there
%are different numbers of time points for each. It reads in two directories
%and overlays the channel with less time points at the matching time and
%each consecutive time frame until another matching time point is found.

%-Jessica Tytell December 2, 2011

%find directories containing dim and bright images
baseFileDir = uigetdir('', 'Select Folder containing images in each time point (base file)');
intermitDir = uigetdir('', 'Select Folder with intermittent images');
outputDir = uigetdir('', 'Select Output Folder for merged images');

%read in files and show in window
baseFiles = imDir(baseFileDir);
disp('baseFiles = ' );
disp(baseFiles);
intermitFiles = imDir(intermitDir);
disp('intermittent Files = ');
disp(intermitFiles);

%write names to own array as string
baseNames = {baseFiles.name};
intermitNames = {intermitFiles.name};

%get length of longer file
nBase = length(baseNames);

%set up incrementer for intermittent files
inc = 1;
for j = 1:nBase
    %get next file and read in
    nextBaseFile = [baseFileDir filesep baseNames{j}];
    baseIm = double(imread(nextBaseFile));
    
    % get time point from filename
    [~, baseBody, baseTime] = getFilenameBody(baseNames{j});
    
    %test if files extensions are the same or smaller
    [~, ~, intTime] = getFilenameBody(intermitNames{inc});
    
    %test to see if next intermittent time point is higher than current
    %base time
    if str2double(baseTime) > str2double(intTime)
        
        inc = inc+1;
        
    end
    
    nextIntFile = [intermitDir filesep intermitNames{inc}];
    intermitIm = double(imread(nextIntFile));
    
    %rgb merge
    rgbIm = ch2rgb(baseIm, intermitIm, []);
    
    outputName= [outputDir filesep baseBody '_merge' baseTime '.TIF'];
    imwrite((rgbIm), outputName, 'tiff');

end

disp('Yay you made it to the end without an error!');
    
    
    
    