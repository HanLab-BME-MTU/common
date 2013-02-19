function  [T_otsu, current_all_matching_bw ]  = geoBasedNmsSeg(nms, classifier_trained)

imageIn=nms;

T_otsu = thresholdOtsu(imageIn);

%%first, get everything
imageMask = imageIn > T_otsu/3;

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

% eliminate the isolated dot and two point lines, just for the ease of
% analysis
% % for i_area = 1 : length(obAreas)
% %     if obAreas(i_area)<=3
% %         labelMask(labelMask == i_area)=0;
% %     end
% % end
%      
% nms_seg_no_brancing = labelMask>0;
%     
% % Label all isolated lines(curves)
% labelMask = bwlabel(nms_seg_no_brancing);
% 
% % Get properties for each of curve
% ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid');
% 
% % Redefine variable for easy of notation
% obAreas = [ob_prop.Area];

nLine = length(obAreas);

obLongaxis = [ob_prop.MajorAxisLength];
obShortaxis = [ob_prop.MinorAxisLength];
obEccentricity = [ob_prop.Eccentricity];

obCentroid = zeros(2, length(obAreas));
obCentroid(:) = [ob_prop.Centroid];

% The ratio of short vs long axis
ratio  = obShortaxis./obLongaxis;

Good_ind=[];
good_lines=cell(1,1);

feature_MinInt = zeros(nLine,1);
feature_MeanInt = zeros(nLine,1);
feature_MaxInt = zeros(nLine,1);

feature_Length = obAreas;
feature_InvCurvature = zeros(nLine,1);

ind_long = find(feature_Length>3);

for i_area = ind_long
%     if obAreas(i_area)>0
        [all_y_i, all_x_i] = find(labelMask == i_area);
        
        INT = imageIn(sub2ind(size(bw_out), round(all_y_i),round(all_x_i)));
%         feature_MinInt(i_area) = min(INT);
        feature_MeanInt(i_area) = mean(INT);
%         feature_MaxInt(i_area) = max(INT);
        
%         bw_i = zeros(size(bw_out));
%         bw_i(sub2ind(size(bw_i), round(all_y_i),round(all_x_i)))=1;
%         end_points_i = bwmorph(bw_i,'endpoints');
%         
%         [y_i, x_i]=find(end_points_i);
%         
%         if isempty(x_i)
%             % if there is no end point, then it is a enclosed circle
%             [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, all_x_i(1),all_y_i(1));
%         else
%             [y_i, x_i]=find(end_points_i);
%             [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
%         end
%         line_smooth_H = fspecial('gaussian',5,1.5);
%         
%         line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
%         line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
%         
%         Vertices = [line_i_x' line_i_y'];
%         Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
%         k=LineCurvature2D(Vertices,Lines);
%         
%         feature_InvCurvature(i_area) = mean(1/k);        
%     end
end


% figure; plot3(feature_Length,feature_MeanInt,feature_InvCurvature,'.');
% 
% feature_Curvature = 1./feature_InvCurvature;
% T_int = thresholdOtsu(imageIn);
% T_L = 6;
% T_K = 0.1;



ind_T_int = 18;
T_xie_int = T_otsu/2*ind_T_int;

Good_ind = find(feature_Length>40);
Bad_ind = [];

for ind_L = 1:40
    T_int = T_otsu/10*ind_T_int;
    
    Good_ind_sub = find(feature_Length'==ind_L & feature_MeanInt>(T_xie_int - ind_L*(T_xie_int/40)));
    Bad_ind_sub = find(feature_Length'==ind_L & feature_MeanInt<=(T_xie_int - ind_L*(T_xie_int/40)));
    
    Good_ind = [Good_ind Good_ind_sub'];
    Bad_ind = [Bad_ind Bad_ind_sub'];
end

[ind_all_y, ind_all_x] = find(labelMask>0);
current_all_matching_bw = zeros(size(labelMask));

for i_E = 1 : length(Good_ind)
    current_good_bw = labelMask==Good_ind(i_E);
    current_all_matching_bw = or(current_all_matching_bw, current_good_bw);
end
