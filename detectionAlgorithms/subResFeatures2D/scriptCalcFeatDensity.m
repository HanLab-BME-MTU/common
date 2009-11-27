
%define where the images are
imageDir = {'/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat02/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat04/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat06/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat08/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat10/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat12/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat14/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat16/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat18/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat20/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat22/',...
    '/home/kj35/.gvfs/orchestra on files.med.harvard.edu/groups/lccb-receptors/Hiro/091019_Q36_SMI_count/Llat24/'};

%define the name of the detection file
detectionFileName = 'detection1.mat';

%define the name of the first image in a series
imageName = '1019_Q36_SMI_count_0001.tif';

%calculate number of movies
numMovies = length(imageDir);

%initialize feature density and its standard deviation
featDensity = NaN(numMovies,2);

%go over all movies
for iMovie = 1 : numMovies
    
    %open the detection file
    load(fullfile(imageDir{iMovie},detectionFileName));
    
    %open the image
    image = imread(fullfile(imageDir{iMovie},imageName));
    
    %find number of nonzero pixels
    cellArea = length(find(image~=0));
    
    %get number of features in each frame
    numFeat = NaN(size(movieInfo));
    for iFrame = 1 : size(movieInfo)
        numFeat(iFrame,1) = size(movieInfo(iFrame).xCoord,1);
    end
    
    %calculate mean and std of feature density
    featDensity(iMovie,:) = [mean(numFeat) std(numFeat)] / cellArea;
    
end