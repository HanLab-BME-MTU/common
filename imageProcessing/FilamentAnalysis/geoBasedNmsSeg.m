function  [T_otsu, current_all_matching_bw ]  = geoBasedNmsSeg(imageIn, classifier_trained, graph_matching_flag)
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

% figure; plot3(feature_Length,feature_MeanInt,feature_InvCurvature,'.');

if(isempty(classifier_trained))
    % if there is no input classifier, build one, with empirical setting
    
    % find the mode of the intensity of the curves/lines
    [hist_n,bin] = hist(feature_MeanInt,200);
    ind_mode = find(hist_n==max(hist_n));
    mode_int = bin(ind_mode(1));
    % And find the Otsu threshold for the intensity
    hotsu = thresholdOtsu(feature_MeanInt(find(feature_MeanInt>mode_int)));
    
    % Set the slanted classification line cutoff as twice of the Otsu with
    % respect to the mode
    T_xie_int =  abs(hotsu - mode_int)*1.5 + mode_int;
        
    % And the length as Otsu threshold
    T_xie_length = 2*max(thresholdOtsu(feature_Length),thresholdRosin(feature_Length));
    
    % Make a classification function as whether it is above the line
    F_classifer = @(i,l) (((T_xie_int + (T_xie_int/T_xie_length)*(-l) )<i));
else
    % when there is an input classifer, use the input one
    F_classifer = classifier_trained;
end

% Select those above the line as our good ones as starting point of graph
% matching
Good_ind = find(F_classifer(feature_MeanInt, feature_Length)>0);

% plot the output image with these good ones
current_all_seg_bw = zeros(size(labelMask));

for i_E = 1 : length(Good_ind)
    current_good_bw = labelMask==Good_ind(i_E);
    current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
end
current_all_matching_bw = current_all_seg_bw;

% graph_matching_flag=1;

% Now if user intended for graph matching part, do it
if(graph_matching_flag==1)
        
    confidency_interval = 0.8;
    current_model = [];
    bw_to_be_matched = current_all_seg_bw;
    
    graph_matching_linking_once; %(current_model, bw_to_be_matched, confidency_interval);
    
    graph_matching_linking_once;
    
    graph_matching_linking_once;
end

