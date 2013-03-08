function  F_classifer  = nms_classifier_train(movieData,imageIn,currentImg)
% nms_classifier_train trains an classifier for segments filaments from input image(nms) based on the geometrical features of the curves/lines in the image and user input
% Input:            
%    MD:                                the MD for the image project
% imageIn:                              the input image, typically the non maximum supress version of the steerable filtering output
% currentImg:                           the greyscale input image

% Output: 
%    F_classifer:                       trained classifier
%

% Liya Ding
% 2013.03



save_tif_flag=1;

if(isempty(imageIn))
    
    % Find the package of Filament Analysis
    nPackage = length(movieData.packages_);
    
    indexFilamentPackage = 0;
    for i = 1 : nPackage
        if(strcmp(movieData.packages_{i}.getName,'FilamentAnalysis')==1)
            indexFilamentPackage = i;
            break;
        end
    end
    
    if(indexFilamentPackage==0)
        msg('Need to be in Filament Package for now.')
        return;
    end
    
    
    nProcesses = length(movieData.processes_);
    
    indexFilamentSegmentationProcess = 0;
    for i = 1 : nProcesses
        if(strcmp(movieData.processes_{i}.getName,'Filament Segmentation')==1)
            indexFilamentSegmentationProcess = i;
            break;
        end
    end
    
    if indexFilamentSegmentationProcess==0
        msg('Please set parameters for steerable filtering.')
        return;
    end
    
    
    funParams=movieData.processes_{indexFilamentSegmentationProcess}.funParams_;
    
    selected_channels = funParams.ChannelIndex;
    StPace_Size = funParams.StPace_Size;
    StPatch_Size = funParams.StPatch_Size;
    Stlowerbound =  funParams.st_lowerbound_localthresholding;
    IntPace_Size = funParams.IntPace_Size;
    IntPatch_Size = funParams.IntPatch_Size;
    Intlowerbound =  funParams.int_lowerbound_localthresholding;
    
    Combine_Way = funParams.Combine_Way;
    Cell_Mask_ind = funParams.Cell_Mask_ind;
    VIF_Outgrowth_Flag = funParams.VIF_Outgrowth_Flag;
    Sub_Sample_Num  = funParams.Sub_Sample_Num;
    
    %% Output Directories
    
    FilamentSegmentationProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'FilamentSegmentation'];
    if (~exist(FilamentSegmentationProcessOutputDir,'dir'))
        mkdir(FilamentSegmentationProcessOutputDir);
    end
    
    
    %%
    indexSteerabeleProcess = 0;
    for i = 1 : nProcesses
        if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
            indexSteerabeleProcess = i;
            break;
        end
    end
    
    if indexSteerabeleProcess==0 && Combine_Way~=2
        msg('Please run steerable filtering first.')
        return;
    end
    
    funParams_st=movieData.processes_{indexSteerabeleProcess}.funParams_;
    
    BaseSteerableFilterSigma = funParams_st.BaseSteerableFilterSigma;
    Levelsofsteerablefilters = funParams_st.Levelsofsteerablefilters;
    ImageFlattenFlag = funParams_st.ImageFlattenFlag;
    
    indexFlattenProcess = 0;
    for i = 1 : nProcesses
        if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
            indexFlattenProcess = i;
            break;
        end
    end
    
    if indexFlattenProcess == 0 && ImageFlattenFlag==2
        display('Please set parameters for Image Flatten.')
        return;
    end
    
    indexCellSegProcess = 0;
    for i = 1 : nProcesses
        if(strcmp(movieData.processes_{i}.getName,'Mask Refinement')==1)
            indexCellSegProcess = i;
            break;
        end
    end
    
    if indexCellSegProcess == 0 && Cell_Mask_ind == 1
        msg('Please run segmentation and refinement first.')
        return;
    end
    
    
    nFrame = movieData.nFrames_;
    
    iChannel = max(selected_channels);
    
    Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  movieData.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel};
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    
    % If steerable filter process is run
    if indexSteerabeleProcess>0
        SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    end
    
    iFrame_index = 1;
    iFrame = 1;
    
    % Read in the intensity image.
    if indexFlattenProcess > 0 && ImageFlattenFlag==2
        currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
    else
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
    end
    currentImg = double(currentImg);
    
    
    load([SteerableChannelOutputDir, filesep, 'steerable_', ...
        filename_short_strs{iFrame},'.mat']);
    
    imageIn = nms;
end

% the threshold defined by Otsu method
T_otsu = thresholdOtsu(imageIn);

% first, get almost all the curves/lines, by using a low threshold
imageMask = imageIn > T_otsu/3;

% further thin it, since the nms version of steerable filtering is not real skeleton
bw_out = bwmorph(imageMask,'thin','inf');

% Find the branching points
nms_seg_brancing = bwmorph(bw_out,'branchpoints');

% Delete these branching points for now
nms_seg_no_brancing = bw_out - nms_seg_brancing;

% again Find the branching points
nms_seg_brancing = bwmorph(nms_seg_no_brancing,'branchpoints');

% Delete these branching points for now
nms_seg_no_brancing = nms_seg_no_brancing - nms_seg_brancing;

% Label all isolated lines(curves)
labelMask = bwlabel(nms_seg_no_brancing);

% Get properties for each of curve
ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid');

% Redefine variable for easy of notation
obAreas = [ob_prop.Area];

nLine = length(obAreas);

% Some feature for later consideration
obLongaxis = [ob_prop.MajorAxisLength];
obShortaxis = [ob_prop.MinorAxisLength];
obEccentricity = [ob_prop.Eccentricity];
obCentroid = zeros(2, length(obAreas));
obCentroid(:) = [ob_prop.Centroid];
% The ratio of short vs long axis
ratio  = obShortaxis./obLongaxis;


feature_MeanInt = nan(nLine,1);
feature_Length = obAreas';

% for the features, only include those curves/lines longer than 4 pixels
ind_long = find(feature_Length>4);

% get the mean intensity of the curves
for i_area = ind_long'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    INT = imageIn(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
    feature_MeanInt(i_area) = mean(INT);
    % this version with the curvature measure, to save time.
end

% figure; plot3(feature_Length,feature_MeanInt,feature_Curvature,'.');


% find the mode of the intensity of the curves/lines
[hist_n,bin] = hist(feature_MeanInt,200);
ind_mode = find(hist_n==max(hist_n));
mode_int = bin(ind_mode(1));
% And find the Otsu threshold for the intensity
hotsu = thresholdOtsu(feature_MeanInt(find(feature_MeanInt>mode_int)));

% the up one
% Set the slanted classification line cutoff as twice of the Otsu with
% respect to the mode
T_xie_int_up =  abs(hotsu - mode_int)*2 + mode_int;

% And the length as Otsu threshold
T_xie_length_up = 2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));

% Make a classification function as whether it is above the line
F_classifer_up = @(i,l) (((T_xie_int_up + (T_xie_int_up/T_xie_length_up)*(-l) )<i));


% the down one
T_xie_int_down =  abs(hotsu - mode_int)*1.5 + mode_int;
T_xie_length_down = 1.5*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));

F_classifer_down = @(i,l) (((T_xie_int_down + (T_xie_int_down/T_xie_length_down)*(-l) )<i));

% the points in between the two classifier is the not sure ones that
% requires annotation
not_sure_ind = find(F_classifer_up(feature_MeanInt, feature_Length)==0 & F_classifer_down(feature_MeanInt, feature_Length)>0);
good_ind = find(F_classifer_up(feature_MeanInt, feature_Length)>0);

F_classifer_notuse =  @(i,l) (((T_xie_int_down/2 + (T_xie_int_down/T_xie_length_down)*(-l) )<i));
bad_ind = find(F_classifer_down(feature_MeanInt, feature_Length)==0 & F_classifer_notuse(feature_MeanInt, feature_Length)>0);
    

h1 = figure(1); hold off;
imagescc(currentImg); hold on;
scrsz = get(0,'ScreenSize');
set(h1,'Position',scrsz);

training_ind = datasample(not_sure_ind,20);

training_bad_ind = [];
training_good_ind = [];

good_bad_label = [];
i_ind = 1;

for i_mark = 1 : 2*length(training_ind)
    i_ind
    iLine = training_ind(i_ind);
    [y,x] = find(labelMask==iLine);
    h1=figure(1);
    plot(x,y,'g.');
    saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_g_',num2str(i_mark),'.jpg']);
    
    AA = [max(1, round(mean(x))-50), min(size(imageIn,2), round(mean(x))+50), ...
        max(1, round(mean(y))-50), min(size(imageIn,1), round(mean(y))+50)];
    axis(AA);
    ch = getkey();
    if(ch==28)
        i_ind = max(1, i_ind - 1);
    else
        i_ind = i_ind + 1;
        if(ch==32)
            good_bad_label(i_ind) = 1;
            plot(x,y,'r.');
        else
            good_bad_label(i_ind) = 2;
            plot(x,y,'b.');
        end
        saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_rb_',num2str(i_mark),'.jpg']);
    
        if(ch==27 || i_ind > length(training_ind))
            i_ind = i_ind-1;
            good_bad_label = good_bad_label(1:length(training_ind));
            break;
        end
        
    end
end

% if everything is good, then lower the cut
if(mean(good_bad_label)>1.8)
    
    % the up one
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int_up =  abs(hotsu - mode_int)*1.5 + mode_int;
    
    % And the length as Otsu threshold
    T_xie_length_up = 1.5*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    % Make a classification function as whether it is above the line
    F_classifer_up = @(i,l) (((T_xie_int_up + (T_xie_int_up/T_xie_length_up)*(-l) )<i));
    
    
    % the down one
    T_xie_int_down =  abs(hotsu - mode_int)*1 + mode_int;
    T_xie_length_down = 1.2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    F_classifer_down = @(i,l) (((T_xie_int_down + (T_xie_int_down/T_xie_length_down)*(-l) )<i));
    
    
    
    
    
    % the points in between the two classifier is the not sure ones that
    % requires annotation
    not_sure_ind = find(F_classifer_up(feature_MeanInt, feature_Length)==0 & F_classifer_down(feature_MeanInt, feature_Length)>0);
    good_ind = find(F_classifer_up(feature_MeanInt, feature_Length)>0);
  
    F_classifer_notuse =  @(i,l) (((T_xie_int_down/2 + (T_xie_int_down/T_xie_length_down)*(-l) )<i));
    bad_ind = find(F_classifer_down(feature_MeanInt, feature_Length)==0 & F_classifer_notuse(feature_MeanInt, feature_Length)>0);
    
    
    h1 = figure(1); hold off;
    imagescc(currentImg); hold on;
    scrsz = get(0,'ScreenSize');
    set(h1,'Position',scrsz);
    
    training_ind = datasample(not_sure_ind,20);
    
    training_bad_ind = [];
    training_good_ind = [];
    
    good_bad_label = [];
    i_ind = 1;
    
    for i_mark = 1 : 2*length(training_ind)
        i_ind
        iLine = training_ind(i_ind);
        [y,x] = find(labelMask==iLine);
        h1=figure(1);
        plot(x,y,'g.');
        saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_g_',num2str(i_mark),'.jpg']);
        
        AA = [max(1, round(mean(x))-50), min(size(imageIn,2), round(mean(x))+50), ...
            max(1, round(mean(y))-50), min(size(imageIn,1), round(mean(y))+50)];
        axis(AA);
        ch = getkey();
        if(ch==28)
            i_ind = max(1, i_ind - 1);
        else
            i_ind = i_ind + 1;
            if(ch==32)
                good_bad_label(i_ind) = 1;
                plot(x,y,'r.');
            else
                good_bad_label(i_ind) = 2;
                plot(x,y,'b.');
            end
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_rb_',num2str(i_mark),'.jpg']);
            
            if(ch==27 || i_ind > length(training_ind))
                i_ind = i_ind-1;
                good_bad_label = good_bad_label(1:length(training_ind));
                break;
            end
            
        end
    end
    
    
end


% if everything is bad, then higher the cut
if(mean(good_bad_label)<1.2)
    
    % the up one
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int_up =  abs(hotsu - mode_int)*3 + mode_int;
    
    % And the length as Otsu threshold
    T_xie_length_up = 3*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    % Make a classification function as whether it is above the line
    F_classifer_up = @(i,l) (((T_xie_int_up + (T_xie_int_up/T_xie_length_up)*(-l) )<i));
    
    
    % the down one
    T_xie_int_down =  abs(hotsu - mode_int)*2 + mode_int;
    T_xie_length_down = 2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    F_classifer_down = @(i,l) (((T_xie_int_down + (T_xie_int_down/T_xie_length_down)*(-l) )<i));
    
    % the points in between the two classifier is the not sure ones that
    % requires annotation
    not_sure_ind = find(F_classifer_up(feature_MeanInt, feature_Length)==0 & F_classifer_down(feature_MeanInt, feature_Length)>0);
    good_ind = find(F_classifer_up(feature_MeanInt, feature_Length)>0);
  
    F_classifer_notuse =  @(i,l) (((T_xie_int_down/2 + (T_xie_int_down/T_xie_length_down)*(-l) )<i));
    bad_ind = find(F_classifer_down(feature_MeanInt, feature_Length)==0 & F_classifer_notuse(feature_MeanInt, feature_Length)>0);
    
    
    h1 = figure(1); hold off;
    imagescc(currentImg); hold on;
    scrsz = get(0,'ScreenSize');
    set(h1,'Position',scrsz);
    
    training_ind = datasample(not_sure_ind,20);
    
    training_bad_ind = [];
    training_good_ind = [];
    
    good_bad_label = [];
    i_ind = 1;
    
    for i_mark = 1 : 2*length(training_ind)
        i_ind
        iLine = training_ind(i_ind);
        [y,x] = find(labelMask==iLine);
        h1=figure(1);
        plot(x,y,'g.');
        saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_g_',num2str(i_mark),'.jpg']);
        
        AA = [max(1, round(mean(x))-50), min(size(imageIn,2), round(mean(x))+50), ...
            max(1, round(mean(y))-50), min(size(imageIn,1), round(mean(y))+50)];
        axis(AA);
        ch = getkey();
        if(ch==28)
            i_ind = max(1, i_ind - 1);
        else
            i_ind = i_ind + 1;
            if(ch==32)
                good_bad_label(i_ind) = 1;
                plot(x,y,'r.');
            else
                good_bad_label(i_ind) = 2;
                plot(x,y,'b.');
            end
            saveas(h1,[FilamentSegmentationChannelOutputDir,'/train_rb_',num2str(i_mark),'.jpg']);
            
            if(ch==27 || i_ind > length(training_ind))
                i_ind = i_ind-1;
                good_bad_label = good_bad_label(1:length(training_ind));
                break;
            end
            
        end
    end
    
    
end


figure(1);
axis auto;

training_good_ind = [good_ind; training_ind(find(good_bad_label(1:end)==1))];
training_bad_ind = [bad_ind; training_ind(find(good_bad_label(1:end)==2))];


train_length_good = feature_Length(training_good_ind);
train_length_bad = feature_Length(training_bad_ind);

train_int_good = feature_MeanInt(training_good_ind);
train_int_bad = feature_MeanInt(training_bad_ind);

% train_cur_good = feature_Curvature(training_good_ind);
% train_cur_bad = feature_Curvature(training_bad_ind);

feature_good = [train_length_good train_int_good ];
feature_bad = [train_length_bad train_int_bad ];

label_good = ones(size(train_length_good));
label_bad = zeros(size(train_length_bad));




feature_training = [feature_Length(training_ind) feature_MeanInt(training_ind)];
label_training = good_bad_label;

mean_good = mean(feature_good);
mean_bad = mean(feature_bad);
v = corr_LDA([feature_good; feature_bad], 2, [length(label_good); length(label_bad)], 0.7);


v*feature_good

v*feature_bad

F_classifer=[];

save('F_classifer.mat','F_classifer');
    


