if retrain_flag==1
    train_mat = [];
    for T_xie_int_grid = T_xie_int_train*(0.8) : (T_xie_int_train*(1.2) - T_xie_int_train*(0.8))/20 : T_xie_int_train*(1.2)
        for T_xie_length_grid = T_xie_length_train*(0.8) : (T_xie_length_train*(1.2) - T_xie_length_train*(0.8))/20 : T_xie_length_train*(1.2)
            
            F_classifer_train = @(i,l) (((T_xie_int_grid + (T_xie_int_grid/T_xie_length_grid)*(-l) -i )));
            train_mat = [train_mat; T_xie_int_grid T_xie_length_grid ...
                (-sum(F_classifer_train(feature_MeanNMS(Matched_ind),...
                feature_Length((Matched_ind))))...
                +sum(F_classifer_train(feature_MeanNMS(UnMatched_ind),...
                feature_Length((UnMatched_ind)))))];
        end
    end
    
    ind = find(train_mat(:,3)==max(train_mat(:,3)));
    
    T_xie_int_train = train_mat(ind(1), 1);
    T_xie_length_train = train_mat(ind(1), 2);
else
    T_xie_int_train = T_xie_int_train*(0.9);
    T_xie_length_train = T_xie_length_train*(0.9);
end

F_classifer = @(int,length) (((T_xie_int_train + (T_xie_int_train/T_xie_length_train)*(-length) )<int));
Good_ind = find(F_classifer(feature_MeanNMS, feature_Length)>0);
Bad_ind = find(F_classifer(feature_MeanNMS, feature_Length)==0);

h12 = figure(12);set(h12,'Visible',set_visible);hold off;
plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');hold on;
plot(feature_Length(Matched_ind),feature_MeanNMS(Matched_ind),'g.');
plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.');


title(['Classifier Plane with matched data before round ',num2str(iIteration)]);
saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_beforeround',num2str(iIteration),'_match_plane.tif']);



% plot the output image with these good ones
current_all_seg_bw = zeros(size(labelMask));
for i_E = 1 : length(Good_ind)
    current_good_bw = labelMask==Good_ind(i_E);
    current_all_seg_bw = or(current_all_seg_bw, current_good_bw);
end

imwrite(double(current_all_seg_bw*3/4),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration+1),'_begin.tif']);


if(exist('h1','var'))
     close(h1);
end
     
[current_model,current_matching_bw, model_ind]...
            = graph_matching_linking_once(current_model, current_all_seg_bw, confidency_interval,imageInt, ...
                                Good_ind,ind_long, model_ind,feature_all,labelMask,master_flassier,iIteration,funParams);
                                
imwrite(double(current_matching_bw/2),[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration+1),'_end.tif']);
h1=figure(1);set(h1,'Visible',set_visible);saveas(h1,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration+1),'_match_color.tif']);
h3=figure(3);set(h3,'Visible',set_visible);saveas(h3,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_round',num2str(iIteration+1),'_all_match_bw.tif']);

good_bw = nms_seg_no_brancing.*current_matching_bw;

original_Matched_ind = Matched_ind;

Matched_ind=[];
UnMatched_ind=[];

for i_area = ind_long'    
    [all_y_i, all_x_i] = find(labelMask == i_area);
    if_matched = mean(good_bw(sub2ind(size(bw_out), round(all_y_i),round(all_x_i))));
    if(if_matched>0.1)
        Matched_ind = [Matched_ind; i_area];
    else
        UnMatched_ind = [UnMatched_ind; i_area];
    end    
end

Good_ind = find(master_flassier(feature_MeanNMS, feature_Length)>0);
Bad_ind = find(master_flassier(feature_MeanNMS, feature_Length)==0);

h12 = figure(12);set(h12,'Visible',set_visible);hold off;
plot(feature_Length(Bad_ind),feature_MeanNMS(Bad_ind),'r.');hold on;
plot(feature_Length(Matched_ind),feature_MeanNMS(Matched_ind),'g.');
plot(feature_Length(setdiff(original_Matched_ind,Good_ind)),feature_MeanNMS(setdiff(original_Matched_ind,Good_ind)),'m*');
plot(feature_Length(Good_ind),feature_MeanNMS(Good_ind),'b.');

title(['Classifier Plane with matched data after round ',num2str(iIteration)]);
saveas(h12,[FilamentSegmentationChannelOutputDir,'/GEO/frame_',num2str(iFrame),'_afterround',num2str(iIteration),'_match_plane.tif']);


% new_curves = zeros(size(imageNMS));
% 
% for ind = (setdiff(Matched_ind,original_Matched_ind))'
% new_curves = imageNMS.*( labelMask==ind);
% title(['ind=',num2str(ind)]);
% % new_curves(1,1)=max(max(imageNMS));
% figure(1);imagescc(new_curves)
% h1 = figure(1);print(h1,'-dtiff',[num2str(ind),'b.tif']);
% 
% end

% new_curves = zeros(size(imageNMS));
% 
% for ind = (setdiff(Matched_ind,original_Matched_ind))'
% new_curves = new_curves+ imageNMS.*( labelMask==ind);
% 
% end
% % new_curves(1,1)=max(max(imageNMS));
% figure(1);imagescc(new_curves)
% h1 = figure(1);print(h1,'-dtiff','brond3_new_all.tif');


