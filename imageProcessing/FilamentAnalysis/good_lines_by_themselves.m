% % function  [level,bw_out] = thresholdLocalSeg(imageIn,choice_of_threshold, level_local_radius, pace, lowerbound, varargin)
load('/home/ld94/files/LCCB/interfil/Liya/fromJess/121127_selected/121127_VimPM60x_D3_02/s2/FilamentAnalysisPackage/FilamentSegmentation/Channel2/DataOutput/steerable_vote_01.mat');
load('/home/ld94/files/LCCB/interfil/Liya/fromJess/121127_selected/121127_VimPM60x_D3_02/s2/FilamentAnalysisPackage/SteerableFiltering/Channel2/steerable_01.mat');
imageIn=nms;

choice_of_threshold = 'Otsu';
level_local_radius = 401;
pace=70;
lowerbound=20;

% level_local_radius = 15;
% pace = 3;
half_pace = round(pace-1)/2;

% Convert to double if necessary
imageIn = double(imageIn);

% Calculate the local and global threshold using the thresholding method
% indicated in choice_of_threshold
[level_img, level_whole] = thresholdLocalCalculate(imageIn,choice_of_threshold,level_local_radius, pace);

% Smooth the threshold map
H = fspecial('gaussian',4*pace+1,pace);
level_img = imfilter(level_img,H,'replicate','same');

% Segmentation using the threshold map
imageMask = imageIn >= level_img/5;

% Get a global mask to eliminate segmentation of detailed noise in the
% background, with a slight lowered threshold to include some boundary
% parts, fill holes and dilate a little to avoid masking off target region

imageMask_whole = imageIn >= level_whole*lowerbound/100;
imageMask_whole = imfill(imageMask_whole,'holes');
% imageMask_whole = imdilate(imageMask_whole,ones(5,5));

% The final segmentation is the intersect of both mask.
bw_out =  imageMask.*imageMask_whole;

% The output level is the gobal threshold 
level = level_whole;
showPlots=1;
if(showPlots==1)
    figure(1); hold off;
    imagesc(imageIn); colormap(gray); hold on
    [y_i, x_i] = find(bw_out > 0);
     plot(x_i,y_i,'.');
%         
%     contour(bw_out,'r')
end

bw_out = bwmorph(bw_out,'thin','inf');

% Find the branching points
nms_seg_brancing = bwmorph(bw_out,'branchpoints');

% Delete these branching points for now
nms_seg_no_brancing = bw_out - nms_seg_brancing;


% Find the branching points
nms_seg_brancing = bwmorph(nms_seg_no_brancing,'branchpoints');

% Delete these branching points for now
nms_seg_no_brancing = nms_seg_no_brancing - nms_seg_brancing;


% Label all isolated lines(curves)
labelMask = bwlabel(nms_seg_no_brancing);

% Get properties for each of curve
ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength','Centroid');

% Redefine variable for easy of notation
obAreas = [ob_prop.Area];
obLongaxis = [ob_prop.MajorAxisLength];
obShortaxis = [ob_prop.MinorAxisLength];
obEccentricity = [ob_prop.Eccentricity];

obCentroid = zeros(2, length(obAreas));
obCentroid(:) = [ob_prop.Centroid];

% The ratio of short vs long axis
ratio  = obShortaxis./obLongaxis;
Good_ind=[];
good_lines=cell(1,1);
for i_area = 1 : length(obAreas)
    if obAreas(i_area)>10
        [y_i, x_i]=find(labelMask == i_area);
        
        bw_i = zeros(size(bw_out));
        bw_i(sub2ind(size(bw_i), round(y_i),round(x_i)))=1;
        
        end_points_i = bwmorph(bw_i,'endpoints');
        
        if isempty(end_points_i)
            % if there is no end point, then it is a enclosed circle
        else
            
            
            [y_i, x_i]=find(end_points_i);
            [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 1000, x_i(1),y_i(1));
               line_smooth_H = fspecial('gaussian',5,1.5);
                
                line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
                line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
             
            Vertices = [line_i_x' line_i_y'];
            Lines=[(1:size(Vertices,1)-1)' (2:size(Vertices,1))'];
            k=LineCurvature2D(Vertices,Lines);
            
            if(mean(abs(k))<0.1)
            
            Good_ind = [Good_ind; i_area];
            
            good_lines{size(Good_ind,1)} = Vertices;
            end
        end
    end
end

figure(1); imagescc(currentImg); hold on;

current_good_bw = zeros(size(bw_out));

current_model = cell(1,1);
count_linking = 0;
current_all_matching_bw = zeros(size(bw_out));

for i_E = 1 : length(Good_ind)
        [y_i, x_i] = find(labelMask == Good_ind(i_E));
        figure(1);
        plot(x_i,y_i,'.');
       
        current_good_bw = labelMask==Good_ind(i_E);
     
        [y_link, x_link] = find(current_good_bw>0);
       
        current_model{i_E} = [x_link, y_link];
       
        current_all_matching_bw = or(current_all_matching_bw, current_good_bw);
       
        figure(2);
        imagescc(current_all_matching_bw);                
     
end


