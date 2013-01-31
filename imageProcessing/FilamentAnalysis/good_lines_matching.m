% Initialize the graph
E = [];
W = [];
angle_array=[];
connection_xy_cell=cell(1,1);

% for now the parameters are set here, without adaptiveness
for i_G = 1 : length(Good_ind)-1
    for j_G = i_G+1 : length(Good_ind)
        
        i_area = Good_ind(i_G);
        j_area = Good_ind(j_G);
        
        try

            x_i =  current_model{i_G}(:,1);
            y_i =  current_model{i_G}(:,2);
            x_j =  current_model{j_G}(:,1);
            y_j =  current_model{j_G}(:,2);
            
            [ind,dist_area(i_area,j_area)] = distance_two_curves([x_i, y_i], [x_j, y_j]);
            
            
            if dist_area(i_area,j_area)<20 && obAreas(i_area)>3 && obAreas(j_area)>3
                
                end_points_i = bwmorph(labelMask == i_area,'endpoints');
                end_points_j = bwmorph(labelMask == j_area,'endpoints');
                
                [y_i, x_i]=find(end_points_i);
                [y_j, x_j]=find(end_points_j);
                [ind, dist_area_endpoint(i_area,j_area)] = distance_two_curves([x_i, y_i], [x_j, y_j]);
                
                [line_i_x, line_i_y] = line_following_with_limit(labelMask == i_area, 10, x_i(ind(1)),y_i(ind(1)));
                [line_j_x, line_j_y] = line_following_with_limit(labelMask == j_area, 10, x_j(ind(2)),y_j(ind(2)));
                
                line_smooth_H = fspecial('gaussian',5,2);
                
                line_i_x = (imfilter(line_i_x, line_smooth_H, 'replicate', 'same'));
                line_i_y = (imfilter(line_i_y, line_smooth_H, 'replicate', 'same'));
                line_j_x = (imfilter(line_j_x, line_smooth_H, 'replicate', 'same'));
                line_j_y = (imfilter(line_j_y, line_smooth_H, 'replicate', 'same'));
                
                
                connect_xy = [line_j_x(1)-line_i_x(1) line_j_y(1)-line_i_y(1)];
                connect_n = connect_xy/norm(connect_xy)/2;
                connect_length = 2*floor(norm(connect_xy))-1;
                
                connect_x = line_i_x(1):connect_n(1):line_i_x(1)+connect_length*connect_n(1);
                connect_y = line_i_y(1):connect_n(2):line_i_y(1)+connect_length*connect_n(2);
                
                
                % flip the i part due to the order is starting the end point
                % toward i, but both connecting and j is leaving i
                
                angle_i = atan2(line_i_x(1)-line_i_x(end),line_i_y(1)-line_i_y(end));
                angle_j = atan2(line_j_x(end)-line_j_x(1),line_j_y(end)-line_j_y(1));
                angle_c = atan2(connect_x(end)-connect_x(1),connect_y(end)-connect_y(1));
                
                
                angle_ij = angle_i - angle_j;
                if angle_ij >pi
                    angle_ij = angle_ij -2*pi;
                end
                if angle_ij < -pi
                    angle_ij = angle_ij + 2*pi;
                end
                if abs(angle_ij) > pi/3
                    continue;
                end
                
                
                angle_ic = angle_i - angle_c;
                if angle_ic >pi
                    angle_ic = angle_ic -2*pi;
                end
                if angle_ic < -pi
                    angle_ic = angle_ic + 2*pi;
                end
                if abs(angle_ic) > pi/3
                    continue;
                end
                
                
                angle_cj = angle_c - angle_j;
                % normalize the angle between -pi to pi
                if angle_cj >pi
                    angle_cj = angle_cj -2*pi;
                end
                if angle_cj < -pi
                    angle_cj = angle_cj + 2*pi;
                end
                if abs(angle_cj) > pi/3
                    continue;
                end
                
                
                
                distance_ij = dist_area_endpoint(i_area,j_area);
                if distance_ij > min(min(obAreas(j_area),obAreas(i_area))/2,20)
                    continue;
                end
                sigma_dis = max(10, min(obAreas(j_area),obAreas(i_area))/4);
                
                sigma_ij = pi/6;
                sigma_ic = pi/6;
                sigma_cj = pi/6;
                
                angle_array = [angle_array; angle_ij angle_ic angle_cj];
                
                W = [W; exp(-(angle_ij.^2)/2/(sigma_ij.^2))...
                    *exp(-(angle_ic.^2)/2/(sigma_ic.^2))...
                    *exp(-(angle_cj.^2)/2/(sigma_cj.^2))...
                    *exp(-(distance_ij.^2)/2/(sigma_dis.^2))...
                    ];
                E = [E; i_area j_area];
                
                connection_xy_cell{size(E,1)} = [connect_x;connect_y];
            end
        end
    end
    
end

M = maxGreedyMatching(E,W);

figure(1); imagescc(currentImg); hold on;

current_matching_bw = zeros(size(bw_out));

current_model = cell(1,1);
count_linking = 0;
current_all_matching_bw = zeros(size(bw_out));

for i_E = 1 : length(W)
    if M(i_E)==1
        
        [y_i, x_i] = find(labelMask == E(i_E,1));
        [y_j, x_j] = find(labelMask == E(i_E,2));
        figure(1);
        plot(x_i,y_i,'.');
        plot(x_j,y_j,'r.');
        
        plot(connection_xy_cell{i_E}(1,:),connection_xy_cell{i_E}(2,:),'y.');
        current_matching_bw = or(labelMask==E(i_E,1),labelMask==E(i_E,2));
        
        current_matching_bw(sub2ind(size(current_matching_bw), round(connection_xy_cell{i_E}(2,:)),round(connection_xy_cell{i_E}(1,:))))=1;
        current_matching_bw = bwmorph(current_matching_bw,'thin','inf');
        
        count_linking = count_linking + 1;
        [y_link, x_link] = find(current_matching_bw>0);
        
        
        current_model{count_linking} = [x_link, y_link];
        
        current_all_matching_bw = or(current_all_matching_bw, current_matching_bw);
        
        figure(2);
        imagescc(current_all_matching_bw);
        
    end
end



for i_E = 1 : length(Good_ind)
        [y_i, x_i] = find(labelMask == Good_ind(i_E));        
        
        % if this good line is not linked to something else
        if(mean(current_all_matching_bw(sub2ind(size(current_matching_bw), round(y_i),round(x_i))))<1)
      
        figure(1); hold on;
        plot(x_i,y_i,'g.');
       
        current_good_bw = labelMask==Good_ind(i_E);
     
        [y_link, x_link] = find(current_good_bw>0);
       
        count_linking = count_linking + 1;       
        current_model{count_linking} = [x_link, y_link];
       
        current_all_matching_bw = or(current_all_matching_bw, current_good_bw);
       
        figure(2);
        imagescc(current_all_matching_bw);
        end
end

