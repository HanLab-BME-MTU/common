function  [T_otsu, current_all_matching_bw, current_model ]  = geoBasedNmsSeg(imageNMS, imageInt, classifier_trained, graph_matching_flag,MaskCell,iFrame)
% geoBasedNmsSeg segments filaments from input image(nms) based on the geometrical features of the curves/lines in the image
% Input:
%    ImageIn:                           the input image, typically the non maximum supress version of the steerable filtering output
%    classifier_trained:                the trained or provided classifier of the curvee, if not provided, use empirical function
% Output:
%    T_otsu:                            the threshold defined by Otsu method for intensity of the input image, just as a format thing here.
%    current_all_matching_bw:           the segmented results, this serves as the starting point of the graphic matching
%

% Liya Ding
% 2013.02

% define a mask to reduce the number of curvelet to classify
if(~isempty(MaskCell) || mean(double(MaskCell))==1)
    
    MaskCell =  imdilate(MaskCell,fspecial('disk', 71)>0);
else
    T_Rosin_otsu = thresholdRosin(imfilter(imageInt,fspecial('gaussian',11,2)));
    MaskCell = imageInt>T_Rosin_otsu*6/3;
    MaskCell = imdilate(MaskCell,fspecial('disk', 71)>0);
end

figure(1); hold off;
subplot(121);
imagescc(imageInt);

% imageNMS = imageNMS.*MaskCell;
% imageInt = imageInt.*MaskCell;
subplot(122);
imagescc(MaskCell);

% pause;
%
% the threshold defined by Otsu method


[hist_n,bin] = hist(imageNMS(find(imageNMS>0)),200);
ind_mode = find(hist_n==max(hist_n));
mode_nms = bin(ind_mode(1));
% And find the Otsu threshold for the intensity
T_otsu = thresholdOtsu(imageNMS(find(imageNMS>mode_nms)));
T_otsu_start =  abs(T_otsu - mode_nms)*0.2 + mode_nms;


imageNMS = imageNMS.*MaskCell;
imageInt = imageInt.*MaskCell;


% first, get almost all the curves/lines, by using a low threshold
imageMask = imageNMS > T_otsu_start;

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
feature_MeanNMS = nan(nLine,1);
feature_Length = obAreas';
smoothed_ordered_points = cell(1,1);
feature_Curvature = nan(nLine,1);

% for the features, only include those curves/lines longer than 4 pixels
ind_long = find(feature_Length>4);

% get the mean intensity of the curves
for i_area = ind_long'
    [all_y_i, all_x_i] = find(labelMask == i_area);
    NMS = imageNMS(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
    feature_MeanNMS(i_area) = mean(NMS);
    INT = imageInt(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
    feature_MeanInt(i_area) = mean(INT);
    % this version with the curvature measure, to save time.
    
    %     bw_i = zeros(size(bw_out));
    %     bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
    %     end_points_i = bwmorph(bw_i,'endpoints');
    %     [y_i, x_i]=find(end_points_i);
    %
    %     if isempty(x_i)
    %         % if there is no end point, then it is a enclosed circle
    %         [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
    %     else
    %         [y_i, x_i]=find(end_points_i);
    %         [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
    %     end
    %
    %     ordered_points{i_area} = [line_i_x, line_i_y];
    %
    %     line_smooth_H = fspecial('gaussian',5,1.5);
    %
    %     line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
    %     line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
    %      smoothed_ordered_points{i_area} = [line_i_x, line_i_y];
    %
    %     Vertices = [line_i_x' line_i_y'];
    %     Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
    %     k=LineCurvature2D(Vertices,Lines);
    %
    %     feature_Curvature(i_area) = mean(k);
end

% figure; plot3(feature_Length,feature_MeanInt,feature_InvCurvature,'.');

if(isempty(classifier_trained))
    % if there is no input classifier, build one, with empirical setting
    
    % find the mode of the intensity of the curves/lines
    [hist_n,bin] = hist(feature_MeanNMS,200);
    ind_mode = find(hist_n==max(hist_n));
    mode_nms = bin(ind_mode(1));
    % And find the Otsu threshold for the intensity
    hotsu = thresholdOtsu(feature_MeanNMS(find(feature_MeanNMS>mode_nms)));
    
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int =  abs(hotsu - mode_nms)*1.0 + mode_nms;
    
    % And the length as Otsu threshold
    T_xie_length = 1.5*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    % Make a classification function as whether it is above the line
    T_xie_int_train = T_xie_int;
    T_xie_length_train = T_xie_length;
    
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    
else
    % when there is an input classifer, use the input one
    F_classifer = classifier_trained;
end

% Select those above the line as our good ones as starting point of graph
% matching
Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);

% plot the output image with these good ones
current_all_seg_bw = zeros(size(labelMask));
current_model = [];

for i_E = 1 : length(Good_ind)
current_good_bw = labelMask==Good_ind(i_E);
current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
[y_i, x_i] = find(labelMask==Good_ind(i_E));

current_model{i_E} = [x_i y_i];
end

current_all_matching_bw = current_all_seg_bw;


% restore parameters from saved trained functions
for int_test = max(feature_MeanNMS):-0.1:0
    if F_classifer(int_test,0)==0;
        T_xie_int_train = int_test;
        break;
    end
end

for length_test = max(feature_Length):-0.1:0
    if F_classifer(0,length_test)==0;
        T_xie_length_train = length_test;
        break;
    end
end

% graph_matching_flag=1;

% Now if user intended for graph matching part, do it
if(graph_matching_flag==1)
    
    confidency_interval = 0.7;
    [current_model,current_matching_bw] = graph_matching_linking_once([], current_all_seg_bw, confidency_interval,imageInt);
      imwrite(double(current_all_seg_bw*3/4),['./GEO/frame_',num2str(iFrame),'_round1_begin.tif']);
   figure(1);close;
    imwrite(double(current_matching_bw/2),['./GEO/frame_',num2str(iFrame),'_round1_end.tif']);
    h1=figure(1);saveas(h1,['./GEO/frame_',num2str(iFrame),'_round1_match_color.tif']);
    
    
    
    T_xie_int_train = T_xie_int_train*0.9;
    T_xie_length_train = T_xie_length_train*0.9;
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
    
    % plot the output image with these good ones
    current_all_seg_bw = zeros(size(labelMask));
    for i_E = 1 : length(Good_ind)
        current_good_bw = labelMask==Good_ind(i_E);
        current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
    end
  imwrite(double(current_all_seg_bw*3/4),['./GEO/frame_',num2str(iFrame),'_round2_begin.tif']);
   figure(1);close;
    [current_model,current_matching_bw] = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt);
    imwrite(double(current_matching_bw/2),['./GEO/frame_',num2str(iFrame),'_round2_end.tif']);
    h1=figure(1);saveas(h1,['./GEO/frame_',num2str(iFrame),'_round2_match_color.tif']);
    
    
    
    T_xie_int_train = T_xie_int_train*0.9;
    T_xie_length_train = T_xie_length_train*0.9;
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
    
    % plot the output image with these good ones
    current_all_seg_bw = zeros(size(labelMask));
    for i_E = 1 : length(Good_ind)
        current_good_bw = labelMask==Good_ind(i_E);
        current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
    end
    
      imwrite(double(current_all_seg_bw*3/4),['./GEO/frame_',num2str(iFrame),'_round3_begin.tif']);
   figure(1);close;
    [current_model,current_matching_bw] = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt);
    imwrite(double(current_matching_bw/2),['./GEO/frame_',num2str(iFrame),'_round3_end.tif']);
    h1=figure(1);saveas(h1,['./GEO/frame_',num2str(iFrame),'_round3_match_color.tif']);
    
    T_xie_int_train = T_xie_int_train*0.9;
    T_xie_length_train = T_xie_length_train*0.9;
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
    
    % plot the output image with these good ones
    current_all_seg_bw = zeros(size(labelMask));
    for i_E = 1 : length(Good_ind)
        current_good_bw = labelMask==Good_ind(i_E);
        current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
    end
  imwrite(double(current_all_seg_bw*3/4),['./GEO/frame_',num2str(iFrame),'_round4_begin.tif']);
   figure(1);close;
    [current_model,current_matching_bw] = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt);
    imwrite(double(current_matching_bw/2),['./GEO/frame_',num2str(iFrame),'_round4_end.tif']);
    h1=figure(1);saveas(h1,['./GEO/frame_',num2str(iFrame),'_round4_match_color.tif']);
    
    T_xie_int_train = T_xie_int_train*0.9;
    T_xie_length_train = T_xie_length_train*0.9;
    F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
    Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
    
    % plot the output image with these good ones
    current_all_seg_bw = zeros(size(labelMask));
    for i_E = 1 : length(Good_ind)
        current_good_bw = labelMask==Good_ind(i_E);
        current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
    end
    imwrite(double(current_all_seg_bw*3/4),['./GEO/frame_',num2str(iFrame),'_round5_begin.tif']);
   figure(1);close;
    [current_model,current_matching_bw] = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt);
    imwrite(double(current_matching_bw/2),['./GEO/frame_',num2str(iFrame),'_round5_end.tif']);
    h1=figure(1);saveas(h1,['./GEO/frame_',num2str(iFrame),'_round5_match_color.tif']);
end

current_all_matching_bw = current_matching_bw>0;