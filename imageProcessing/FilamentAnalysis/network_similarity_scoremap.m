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
angle_map_1_2_A = nan(img_size);
angle_map_2_1_A = nan(img_size);
angle_map_1_2_B = nan(img_size);
angle_map_2_1_B = nan(img_size);

  [MT_digital_model,MT_orientation_model,MT_XX,MT_YY,MT_OO] ...
        = filament_model_to_digital_with_orientation(MT_current_model);
  [VIF_digital_model,VIF_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
        = filament_model_to_digital_with_orientation(VIF_current_model);

[Y1, X1] = find(VIF_current_seg>0);
[Y2, X2] = find(MT_current_seg>0);


X1=VIF_XX;
Y1=VIF_YY;

X2=MT_XX;
Y2=MT_YY;

O1 = VIF_OO;
O2 = MT_OO;


% find in the pool of the points in the 2nd channel, the closest point of
% the points in the 1st channel, so query is 1, pool is 2.
[idx_cell, dist_cell] = KDTreeBallQuery([Y2 X2],[Y1 X1],radius);

for iQ = 1 : length(Y1)
    idx_this = idx_cell{iQ};
    dist_this = dist_cell{iQ};
    
    [sort_dist,sort_IX] = sort(dist_this);
    sort_idx = idx_this(sort_IX);
    
    distance_map_1_2(sub2ind(img_size,Y1(iQ),X1(iQ)))=sort_dist(1);
    angle_map_1_2_A(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
        VIF_OO(iQ) - MT_OO(sort_idx(1));
    if sort_dist(1)==sort_dist(2) && Y2(sort_idx(1))==Y2(sort_idx(2))...
            && X2(sort_idx(1))==X2(sort_idx(2))
        angle_map_1_2_B(sub2ind(img_size,Y1(iQ),X1(iQ))) = ...
            VIF_OO(iQ) - MT_OO(sort_idx(2));
    end
end

angle_map_1_2_A(angle_map_1_2_A>pi/2) = angle_map_1_2_A(angle_map_1_2_A>pi/2) -pi;
angle_map_1_2_A(angle_map_1_2_A<-pi/2) = angle_map_1_2_A(angle_map_1_2_A<-pi/2) +pi;

angle_map_1_2_A(angle_map_1_2_A>pi/2) = angle_map_1_2_A(angle_map_1_2_A>pi/2) -pi;
angle_map_1_2_A(angle_map_1_2_A<-pi/2) = angle_map_1_2_A(angle_map_1_2_A<-pi/2) +pi;

angle_map_1_2_B(angle_map_1_2_B>pi/2) = angle_map_1_2_B(angle_map_1_2_B>pi/2) -pi;
angle_map_1_2_B(angle_map_1_2_B<-pi/2) = angle_map_1_2_B(angle_map_1_2_B<-pi/2) +pi;

angle_map_1_2_B(angle_map_1_2_B>pi/2) = angle_map_1_2_B(angle_map_1_2_B>pi/2) -pi;
angle_map_1_2_B(angle_map_1_2_B<-pi/2) = angle_map_1_2_B(angle_map_1_2_B<-pi/2) +pi;

angle_map_1_2 = min(abs(angle_map_1_2_A),abs(angle_map_1_2_B));

% % if the dist is farther than radius, gate it, so this index is for 1
% ind_gate = find(dist<=radius);
% 
% % then update the index for founded closest points
% idx_close = idx(ind_gate);
% 
% % screen out the gated ones
% X1_close = X1(ind_gate);
% Y1_close = Y1(ind_gate);
% X2_close = X2(idx_close);
% Y2_close = Y2(idx_close);
% 
% dist_gate = dist(ind_gate);
% 
% distance_map_1_2(sub2ind(img_size,Y1_close,X1_close))=dist_gate;
% angle_map_1_2(sub2ind(img_size,Y1_close,X1_close)) = ...
%     VIF_orientation(sub2ind(img_size,Y1_close,X1_close)) ...
%     - MT_orientation(sub2ind(img_size,Y2_close,X2_close));
% 

 %%
%then the other way around 2->1


[idx_cell, dist_cell] = KDTreeBallQuery([Y1 X1],[Y2 X2],radius);

for iQ = 1 : length(Y2)
    idx_this = idx_cell{iQ};
    dist_this = dist_cell{iQ};
    
    [sort_dist,sort_IX] = sort(dist_this);
    sort_idx = idx_this(sort_IX);
    if(length(sort_dist)>0)
        distance_map_2_1(sub2ind(img_size,Y2(iQ),X2(iQ)))=sort_dist(1);
        angle_map_2_1_A(sub2ind(img_size,Y2(iQ),X2(iQ))) = ...
            MT_OO(iQ) - VIF_OO(sort_idx(1));
    end
    
    if(length(sort_dist)>1)        
        if sort_dist(1)==sort_dist(2) && Y1(sort_idx(1))==Y1(sort_idx(2))...
                && X1(sort_idx(1))==X1(sort_idx(2))
            angle_map_2_1_B(sub2ind(img_size,Y2(iQ),X2(iQ))) = ...
                MT_OO(iQ) - VIF_OO(sort_idx(2));
        end
    end
end

angle_map_2_1_A(angle_map_2_1_A>pi/2) = angle_map_2_1_A(angle_map_2_1_A>pi/2) -pi;
angle_map_2_1_A(angle_map_2_1_A<-pi/2) = angle_map_2_1_A(angle_map_2_1_A<-pi/2) +pi;

angle_map_2_1_A(angle_map_2_1_A>pi/2) = angle_map_2_1_A(angle_map_2_1_A>pi/2) -pi;
angle_map_2_1_A(angle_map_2_1_A<-pi/2) = angle_map_2_1_A(angle_map_2_1_A<-pi/2) +pi;

angle_map_2_1_B(angle_map_2_1_B>pi/2) = angle_map_2_1_B(angle_map_2_1_B>pi/2) -pi;
angle_map_2_1_B(angle_map_2_1_B<-pi/2) = angle_map_2_1_B(angle_map_2_1_B<-pi/2) +pi;

angle_map_2_1_B(angle_map_2_1_B>pi/2) = angle_map_2_1_B(angle_map_2_1_B>pi/2) -pi;
angle_map_2_1_B(angle_map_2_1_B<-pi/2) = angle_map_2_1_B(angle_map_2_1_B<-pi/2) +pi;

angle_map_2_1 = min(abs(angle_map_2_1_A),abs(angle_map_2_1_B));
% 
% 
% % find in the pool of the points in the 1nd channel, the closest point of
% % the points in the 2st channel, so query is 2, pool is 1 .
% [idx, dist] = KDTreeClosestPoint([Y1 X1],[Y2 X2]);
% 
% % if the dist is farther than radius, gate it, so this index is for 2
% ind_gate = find(dist<=radius);
% 
% % then update the index for founded closest points
% idx_close = idx(ind_gate);
% 
% % screen out the gated ones
% X2_close = X2(ind_gate);
% Y2_close = Y2(ind_gate);
% X1_close = X1(idx_close);
% Y1_close = Y1(idx_close);
% 
% dist_gate = dist(ind_gate);
% 
% distance_map_2_1(sub2ind(img_size,Y2_close,X2_close)) = dist_gate;
% angle_map_2_1(sub2ind(img_size,Y2_close,X2_close)) = ...
%     MT_orientation(sub2ind(img_size,Y2_close,X2_close)) ...
%     - VIF_orientation(sub2ind(img_size,Y1_close,X1_close));
% 
% 
% angle_map_1_2(angle_map_1_2>pi/2) = angle_map_1_2(angle_map_1_2>pi/2) -pi;
% angle_map_1_2(angle_map_1_2<-pi/2) = angle_map_1_2(angle_map_1_2<-pi/2) +pi;
% 
% 
% angle_map_2_1(angle_map_2_1>pi/2) = angle_map_2_1(angle_map_2_1>pi/2) -pi;
% angle_map_2_1(angle_map_2_1<-pi/2) = angle_map_2_1(angle_map_2_1<-pi/2) +pi;
% 
% 
% 
% imwrite(distance_map_1_2,[outdir,filesep,'Dis12_frame_',num2str(iFrame),'.tif']);
% imwrite(distance_map_2_1,[outdir,filesep,'Dis21_frame_',num2str(iFrame),'.tif']);
% imwrite(angle_map_2_1,[outdir,filesep,'Ang21_frame_',num2str(iFrame),'.tif']);
% imwrite(angle_map_1_2,[outdir,filesep,'Ang12_frame_',num2str(iFrame),'.tif']);
% 

h6=figure(6);
imagesc_nan_neg(distance_map_1_2,radius*1);axis equal;axis off;
axis([150 250 140 200]);saveas(h6,[outdir,filesep,'Dis12_frame_',num2str(iFrame),'.tif']);

h6=figure(6);
imagesc_nan_neg(distance_map_2_1,radius*1);axis equal;axis off;
axis([150 250 140 200]);saveas(h6,[outdir,filesep,'Dis21_frame_',num2str(iFrame),'.tif']);

h6=figure(6);
imagesc_nan_neg(angle_map_2_1,-pi/(3));axis equal;axis off;
axis([150 250 140 200]);saveas(h6,[outdir,filesep,'Ang12_frame_',num2str(iFrame),'.tif']);

h6=figure(6);
imagesc_nan_neg(angle_map_1_2,-pi/(3));axis equal;axis off;
axis([150 250 140 200]);
saveas(h6,[outdir,filesep,'Ang21_frame_',num2str(iFrame),'.tif']);

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
Weight_mask = fspecial('gaussian',2*radius+1,radius/(2));
Weight_mask = Weight_mask(sub2ind([2*radius+1,2*radius+1,],cy,cx));

% for all these points
for j = 1 : length(Y)
%     j
    x = X(j);
    y = Y(j);
        
    dis_1 = distance_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));
    dis_2 = distance_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1));

    ang_1 = abs(angle_map_1_2_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1)));
    ang_2 = abs(angle_map_2_1_pad(sub2ind(img_size+2*radius,y+cy-1,x+cx-1)));
    
    if(~isempty(dis_1))
        score_maps_distance_1_2(y,x) = sum(dis_1(~isnan(dis_1)).*Weight_mask(~isnan(dis_1)))/sum(Weight_mask(~isnan(dis_1)));
        score_maps_angle_1_2(y,x) = sum(ang_1(~isnan(ang_1)).*Weight_mask(~isnan(ang_1)))/sum(Weight_mask(~isnan(ang_1)));
    end
    if(~isempty(dis_2))
        score_maps_distance_2_1(y,x) = sum(dis_2(~isnan(dis_2)).*Weight_mask(~isnan(dis_2)))/sum(Weight_mask(~isnan(dis_2)));       
        score_maps_angle_2_1(y,x) = sum(ang_2(~isnan(ang_2)).*Weight_mask(~isnan(ang_2)))/sum(Weight_mask(~isnan(ang_2)));
    end    
end

score_maps_distance_2_1(isnan(score_maps_distance_2_1))=radius;
score_maps_distance_1_2(isnan(score_maps_distance_1_2))=radius;
score_maps_angle_2_1(isnan(score_maps_angle_2_1))=pi/2;
score_maps_angle_1_2(isnan(score_maps_angle_1_2))=pi/2;

h6=figure(6);
imagesc_nan_neg(score_maps_distance_1_2,0);axis equal;axis off;
saveas(h6,[outdir,filesep,'LVDis12_frame_',num2str(iFrame),'.tif']);

h6=figure(6);
imagesc_nan_neg(score_maps_distance_2_1,0);axis equal;axis off;
saveas(h6,[outdir,filesep,'LVDis21_frame_',num2str(iFrame),'.tif']);

h6=figure(6);
imagesc_nan_neg(score_maps_angle_1_2,0);axis equal;axis off;
saveas(h6,[outdir,filesep,'LVAng12_frame_',num2str(iFrame),'.tif']);

h6=figure(6);
imagesc_nan_neg(score_maps_angle_2_1,0);axis equal;axis off;
saveas(h6,[outdir,filesep,'LVAng21_frame_',num2str(iFrame),'.tif']);



h3=figure(3); imagesc_nan_neg(score_maps_distance_2_1+score_maps_distance_1_2,0);axis equal;axis off;
saveas(h3,[outdir,filesep,'VIFMT_dis_frame_',num2str(iFrame),'.tif']);

saveas(h3,[outdir,filesep,'VIFMT_dis_frame_',num2str(iFrame),'.fig']);

h4=figure(4); imagesc_nan_neg(abs(score_maps_angle_2_1/2)+abs(score_maps_angle_1_2/2),0);axis equal;axis off;
saveas(h4,[outdir,filesep,'VIFMT_ang_frame_',num2str(iFrame),'.tif']);

saveas(h4,[outdir,filesep,'VIFMT_ang_frame_',num2str(iFrame),'.fig']);

similarity_score = exp(-(score_maps_distance_2_1+score_maps_distance_1_2).^2/((radius*1.5)^2))...
    .*exp(-(abs(score_maps_angle_2_1/2)+abs(score_maps_angle_1_2/2)).^2/((pi/4)^2));
save([outdir,filesep,'VIFMT_sm_maps_frame_',num2str(iFrame),'.mat'], ...
    'score_maps_distance_2_1','score_maps_distance_1_2',...
    'score_maps_angle_2_1','score_maps_angle_1_2',...
    'similarity_score');

similarity_score(similarity_score<0.2)=0.2;
h5=figure(5); imagesc_nan_neg(similarity_score,0);axis equal;axis off;
saveas(h5,[outdir,filesep,'VIFMT_sm_score_frame_',num2str(iFrame),'.tif']);

saveas(h5,[outdir,filesep,'VIFMT_sm_score_frame_',num2str(iFrame),'.fig']);


%%
% Crossing comparison

[idx_cell, dist_cell] = KDTreeBallQuery([Y1 X1],[Y1 X1],0.1);

crossing_flag_1 = zeros(1, length(X1));
angle_crossing_1 = nan(1, length(X1));

for iQ = 1 : length(Y1)
    idx_this = idx_cell{iQ};
    if(length(idx_this)>1)
        
        AA = O1(idx_this(1))-O1(idx_this(2));
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        
        if(abs(AA)>0.2)
            crossing_flag_1(iQ)=1;
            angle_crossing_1(iQ)=O1(idx_this(1))-O1(idx_this(2));
            ind_pair(iQ,1)=idx_this(1);
            ind_pair(iQ,2)=idx_this(2);
        end
    end
end
angle_crossing_1(angle_crossing_1>pi/2) = angle_crossing_1(angle_crossing_1>pi/2)-pi;
angle_crossing_1(angle_crossing_1<-pi/2) = angle_crossing_1(angle_crossing_1<-pi/2)+pi;
angle_crossing_1(angle_crossing_1>pi/2) = angle_crossing_1(angle_crossing_1>pi/2)-pi;
angle_crossing_1(angle_crossing_1<-pi/2) = angle_crossing_1(angle_crossing_1<-pi/2)+pi;



[idx_cell, dist_cell] = KDTreeBallQuery([Y2 X2],[Y2 X2],0.2);

crossing_flag_2 = zeros(1, length(X2));
angle_crossing_2 = nan(1, length(X2));

for iQ = 1 : length(Y2)
    idx_this = idx_cell{iQ};
    if(length(idx_this)>1)
        
        
         AA = O2(idx_this(1))-O2(idx_this(2));
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        
        if(abs(AA)>0.2)
            crossing_flag_2(iQ)=1;
            angle_crossing_2(iQ)=O2(idx_this(2))-O2(idx_this(2));     
        end
    end    
end
angle_crossing_2(angle_crossing_2>pi/2) = angle_crossing_2(angle_crossing_2>pi/2)-pi;
angle_crossing_2(angle_crossing_2<-pi/2) = angle_crossing_2(angle_crossing_2<-pi/2)+pi;
angle_crossing_2(angle_crossing_2>pi/2) = angle_crossing_2(angle_crossing_2>pi/2)-pi;
angle_crossing_2(angle_crossing_2<-pi/2) = angle_crossing_2(angle_crossing_2<-pi/2)+pi;



[idx, dist] = KDTreeClosestPoint([Y2(crossing_flag_2>0) X2(crossing_flag_2>0)],...
    [Y1(crossing_flag_1>0) X1(crossing_flag_1>0) ]);

dist_pool_for_crossing = [dist_pool_for_crossing;dist];
ang_pool_for_crossing = [ang_pool_for_crossing;(angle_crossing_1(crossing_flag_1>0))'];


[idx, dist] = KDTreeClosestPoint([Y1(crossing_flag_1>0) X1(crossing_flag_1>0)],...
[Y2(crossing_flag_2>0) X2(crossing_flag_2>0)]);

dist_pool_for_crossing = [dist_pool_for_crossing;dist];
ang_pool_for_crossing = [ang_pool_for_crossing;(angle_crossing_2(crossing_flag_2>0))'];
