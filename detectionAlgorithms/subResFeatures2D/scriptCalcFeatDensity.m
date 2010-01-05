
%define where the images are
imageDir = {'/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_02/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_04/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_06/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_08/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_10/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_clnb_12/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_bleb_02/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_bleb_05/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_bleb_07/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_bleb_09/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_noco_02/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_noco_04/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_noco_06/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_noco_11/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_noco_17/images/',...
    '/orchestra/groups/lccb-receptors/Hiro/071229_many/control/con_noco_21/images/'};

%define the name of the detection file
detectionFileName = 'detection1.mat';

%define the name of the first image in a series
imageName = '1126Fab36_SPT_jasp_0001.tif';

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