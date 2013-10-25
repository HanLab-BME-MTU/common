function   current_seg = filament_model_to_seg_bwim(current_model,img_size,break_list)

model_length = length(current_model);

current_seg =zeros(img_size);

for iFila = 1 : model_length    
    line_i_x = current_model{iFila}(:,1);
    line_i_y = current_model{iFila}(:,2);
    
    current_seg(sub2ind(img_size, round(line_i_y),round(line_i_x)))=1;    
    
end
  