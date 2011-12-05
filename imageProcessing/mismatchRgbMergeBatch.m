function mismatchRgbMergeBatch(baseFileDir, intermitDir)
%This function is designed to make RGB images of two channels where there
%are different numbers of time points for each. It reads in two directories
%and overlays the channel with less time points at the matching time and
%each consecutive time frame until another matching time point is found.

%-Jessica Tytell December 2, 2011

% parse input
ip = inputParser;
ip.addRequired('baseFileDir', @ischar); 
ip.addRequired('intermitDir', @ischar);

%set location and make output dir
outputLocation = fileparts(dimDir);
outputDir = ([outputLocation filesep 'VimentinPlusTipMerge']);
mkdir(outputDir);

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

%add waitbar for impatient people (or people who don't trust this code)
h = waitbar(0,'Please wait...');

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
    
    %update waitbar
    waitbar(j / nBase)

end
close(h);

disp('<snoopy dance>');
    
    
    
    