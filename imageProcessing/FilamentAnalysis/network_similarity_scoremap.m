% function score_maps = network_similarity_scoremap(VIF_orientation, VIF_current_model,VIF_current_seg, ...
%     MT_orientation,MT_current_model,MT_current_seg,radius)
% function for calculation the similarity of two networks
% Liya Ding 06.2013.
radius=20;

score_maps = cell(1,7);
img_size = size(VIF_current_seg);

distance_map_1_2 = nan(img_size);
distance_map_2_1 = nan(img_size);
angle_map_1_2 = nan(img_size);
angle_map_2_1 = nan(img_size);

[Y1, X1] = find(VIF_current_seg>0);
[Y2, X2] = find(MT_current_seg>0);

% find in the pool of the points in the 2nd channel, the closest point of
% the points in the 1st channel, so query is 1, pool is 2.
[idx, dist] = KDTreeClosestPoint([Y2 X2],[Y1 X1]);

% if the dist is farther than radius, gate it, so this index is for 1
ind_gate = find(dist<=radius);

% then update the index for founded closest points
idx_close = idx(ind_gate);

% screen out the gated ones
X1_close = X1(ind_gate);
Y1_close = Y1(ind_gate);
X2_close = X2(idx_close);
Y2_close = Y2(idx_close);

dist_gate = dist(ind_gate);

distance_map_1_2(sub2ind(img_size,Y1_close,X1_close))=dist_gate;
angle_map_1_2(sub2ind(img_size,Y1_close,X1_close)) = ...
    VIF_orientation(sub2ind(img_size,Y1_close,X1_close)) ...
    - MT_orientation(sub2ind(img_size,Y2_close,X2_close));

% then the other way around 2->1



% find in the pool of the points in the 1nd channel, the closest point of
% the points in the 2st channel, so query is 2, pool is 1 .
[idx, dist] = KDTreeClosestPoint([Y1 X1],[Y2 X2]);

% if the dist is farther than radius, gate it, so this index is for 2
ind_gate = find(dist<=radius);

% then update the index for founded closest points
idx_close = idx(ind_gate);

% screen out the gated ones
X2_close = X2(ind_gate);
Y2_close = Y2(ind_gate);
X1_close = X1(idx_close);
Y1_close = Y1(idx_close);

dist_gate = dist(ind_gate);

distance_map_2_1(sub2ind(img_size,Y2_close,X2_close)) = dist_gate;
angle_map_2_1(sub2ind(img_size,Y2_close,X2_close)) = ...
    MT_orientation(sub2ind(img_size,Y2_close,X2_close)) ...
    - VIF_orientation(sub2ind(img_size,Y1_close,X1_close));


angle_map_1_2(angle_map_1_2>pi/2) = angle_map_1_2(angle_map_1_2>pi/2) -pi;
angle_map_1_2(angle_map_1_2<-pi/2) = angle_map_1_2(angle_map_1_2<-pi/2) +pi;


angle_map_2_1(angle_map_2_1>pi/2) = angle_map_2_1(angle_map_2_1>pi/2) -pi;
angle_map_2_1(angle_map_2_1<-pi/2) = angle_map_2_1(angle_map_2_1<-pi/2) +pi;

whole_ROI = imdilate(VIF_current_seg,ones(radius,radius)) + imdilate(MT_current_seg,ones(radius,radius))>0;
% 
% % disgard the boundary ones
% whole_ROI(1:radius+1,:)=0;
% whole_ROI(end-radius:end,:)=0;
% whole_ROI(:,1:radius+1)=0;
% whole_ROI(:,end-radius:end)=0;

% find the points of interest
[Y,X] = find(whole_ROI>0);

score_maps_distance_1_2 = nan(img_size);
score_maps_distance_2_1 = nan(img_size);
score_maps_angle_1_2 = nan(img_size);
score_maps_angle_2_1 = nan(img_size);

[cy,cx] = find(fspecial('disk',radius)>0);


distance_map_1_2_pad = nan(img_size+2*radius);
distance_map_2_1_pad = nan(img_size+2*radius);
angle_map_1_2_pad = nan(img_size+2*radius);
angle_map_2_1_pad = nan(img_size+2*radius);

distance_map_1_2_pad(radius+1:end-radius,radius+1:end-radius) = distance_map_1_2;
distance_map_2_1_pad(radius+1:end-radius,radius+1:end-radius) = distance_map_2_1;
angle_map_1_2_pad(radius+1:end-radius,radius+1:end-radius) = angle_map_1_2;
angle_map_2_1_pad(radius+1:end-radius,radius+1:end-radius) = angle_map_2_1;

% 
% distance_map_1_2_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% distance_map_2_1_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% angle_map_1_2_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% angle_map_2_1_pad_stack = nan(img_size(1)+2*radius,img_size(2)+2*radius,length(cx));
% 
% 
% for s = 1 : length(cx)
%     distance_map_1_2_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = distance_map_1_2;
%     distance_map_2_1_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = distance_map_2_1;
%     angle_map_1_2_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = angle_map_1_2;
%     angle_map_2_1_pad_stack(cy(s):cy(s)+img_size(1)-1,cx(s):cx(s)+img_size(2)-1, s) = angle_map_2_1;
% end
% 
% 


% for all these points
for j = 1 : length(Y)
    j
    x = X(j);
    y = Y(j);
        
    dis_1 = distance_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));
    dis_2 = distance_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));

    ang_1 = angle_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));
    ang_2 = angle_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));
    
    if(~isempty(dis_1))
        score_maps_distance_1_2(y,x) = mean(dis_1(~isnan(dis_1)));
        score_maps_angle_1_2(y,x) = mean(ang_1(~isnan(ang_1)));
    end
    if(~isempty(dis_2))
        score_maps_distance_2_1(y,x) = mean(dis_2(~isnan(dis_2)));        
        score_maps_angle_2_1(y,x) = mean(ang_2(~isnan(ang_2)));
    end    
end

score_maps_distance_2_1(isnan(score_maps_distance_2_1))=radius;
score_maps_distance_1_2(isnan(score_maps_distance_1_2))=radius;
score_maps_angle_2_1(isnan(score_maps_angle_2_1))=pi/2;
score_maps_angle_1_2(isnan(score_maps_angle_1_2))=pi/2;

h3=figure(3); imagesc(score_maps_distance_2_1+score_maps_distance_1_2);axis equal;axis off;
saveas(h3,[outdir,filesep,'MTMT_dis_frame_',num2str(iFrame),'.tif']);

h4=figure(4); imagesc(abs(score_maps_angle_2_1/2)+abs(score_maps_angle_1_2/2));axis equal;axis off;
saveas(h4,[outdir,filesep,'MTMT_ang_frame_',num2str(iFrame),'.tif']);

similarity_score = exp(-(score_maps_distance_2_1+score_maps_distance_1_2).^2/((radius)^2))...
    .*exp(-(abs(score_maps_angle_2_1/2)+abs(score_maps_angle_1_2/2)).^2/((pi/6)^2));
save([outdir,filesep,'MTMT_sm_maps_frame_',num2str(iFrame),'.mat'], ...
    'score_maps_distance_2_1','score_maps_distance_1_2',...
    'score_maps_angle_2_1','score_maps_angle_1_2',...
    'similarity_score');

similarity_score(similarity_score<0.3)=0.3;
h5=figure(5); imagesc(similarity_score);axis equal;axis off;
saveas(h5,[outdir,filesep,'MTMT_sm_score_frame_',num2str(iFrame),'.tif']);



