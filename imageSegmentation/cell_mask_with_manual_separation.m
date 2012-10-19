function cell_mask_with_manual_separation(MD,startframe)
% manually help the cell segmentation of a single cell out off many
% touching or not touching cells

% Input: MD: the movieData object

% Output: nothing as return value, but the segmented images are save to
% refinement results folder under folder name "Cell_Single"

% Usage: for each frame, two step
% 1. Click at the cell of interest, or starting the second frame, just keep
% checking the cell at the same location as previous frame
% 2. Cut off any connecting part, can go as many times as needed. User
% input the line of cutting by clicking the points along it, the program is
% using the poly ROI definition input, NOTE here the line connecting the
% first and last points will be ignored(since we are looking for a cutting line),
% so you don't need to worry about circling the whole cell into the ROI.
% When finish or you think there is nothing to cut of, double click as if the 
% ROI is just one point, the code will detect you finish cutting.

% Liya Ding 2012.10.
% low budget code for helping the cell segmentation, well, cell separation.

%% Find the process of segmentation and mask refinement.
nProcesses = length(MD.processes_);

indexMBSegProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Thresholding')==1)
        indexMBSegProcess = i;
        break;
    end
end

if indexMBSegProcess==0
    msg('Please run segmentation and refinement first.')
    return;
end

indexMBSegRefineProcess = 0;
for i = 1 : nProcesses
    if(strcmp(MD.processes_{i}.getName,'Mask Refinement')==1)
        indexMBSegRefineProcess = i;
        break;
    end
end

if indexMBSegRefineProcess==0
    msg('Please run segmentation and refinement first.')
    return;
end


indexMBChannel = 1;
nFrame = MD.nFrames_;

%% New output folder
Cell_manual_seg_folder = [MD.processes_{indexMBSegRefineProcess}.outFilePaths_{indexMBChannel},filesep,'Cell_Single'];
mkdir(Cell_manual_seg_folder);

%%
H_close = fspecial('disk',2);
H_close = H_close>0;

%% Control arrays
poisition_flag=zeros(1,nFrame);
manual_finish=zeros(1,nFrame);

%% For each frame, do manual helping/checking
for iFrame = startframe : nFrame
    % load segmented mask ( not refined ones)
    MaskMB = MD.processes_{indexMBSegProcess}.loadChannelOutput(indexMBChannel,iFrame);
    % load images, boost them for visualization purpose
    currentImg = MD.channels_(indexMBChannel).loadImage(iFrame);
    currentImg = sqrt(currentImg/max(max(currentImg)));
    currentImg = sqrt(currentImg/max(max(currentImg)));
    
    % For display
    MaskMB_display = currentImg*250;
    MaskMB_display(MaskMB>0)=255;
      
    %% Step 1: select cell
    
    % Show the mask and the boundary
    figure(1); hold off;
    imagesc(MaskMB_display);  axis off; axis image; colormap(gray);
    RoiYX = bwboundaries(MaskMB);
    contourYX = RoiYX{1};
    hold on;
    plot(contourYX(:,2),contourYX(:,1),'r');
    
    if iFrame == startframe
        height_image = size(MaskMB,1);
        width_image = size(MaskMB,2);
       
        title(['Frame: ',num2str(iFrame),', please click on which cell you want to help segment']);
        % For the first frame, the user has to define a valid cell center
        while poisition_flag(iFrame) ==0
            [x,y,button]  = ginput(1);
            if(button==1)
                if(x > 0 && x <=width_image && y >0 && y <= height_image)
                    poisition_flag(iFrame) = 1;
                    x= round(x); y=round(y);
                    user_center_y(iFrame)=y;
                    user_center_x(iFrame)=x;
                end
            end
        end
    else
        title(['Frame: ',num2str(iFrame),', if cell is still where the dot is, enter or right lick. Otherwise, click on the cell.']);
        % For the later frame, user can define cell center, or just keep
        % the previous one
        while poisition_flag(iFrame) ==0
            [x,y,button]  = ginput(1);
            if(button==1)
                if(x > 0 && x <=width_image && y >0 && y <= height_image)
                    x= round(x); y=round(y);
                    starting_poisition_flag = 1;
                    user_center_y(iFrame)=y;
                    user_center_x(iFrame)=x;
                    poisition_flag(iFrame)=1;
                end
            else
                user_center_y(iFrame)=user_center_y(iFrame-1);
                user_center_x(iFrame)=user_center_x(iFrame-1);
                poisition_flag(iFrame)=1;
            end
        end
    end
    
    % Get rid of the cells not touching to the current selected cell first
    labelMask = bwlabel(MaskMB,4);
    MaskMB_new = labelMask==labelMask(round(y),round(x));
    user_center_y(iFrame)=y;
    user_center_x(iFrame)=x;
    
    MaskMB_new = imerode(MaskMB_new,H_close);
    MaskMB_new = imdilate(MaskMB_new,H_close);
    
    % Keep largest piece
    labelMask = bwlabel(MaskMB_new,4);
    obAreas = regionprops(labelMask,'Area');
    
    if length(obAreas) > 1
        obAreas = [obAreas.Area];
        [dummy,iSort] = sort(obAreas,'descend');
        MaskMB_new = labelMask == iSort(1);
    end
    % Fill the hole
    MaskMB_new = imfill(MaskMB_new,'hole');
    
       
    %% Step 2: Cutting touching parts
    
    while(manual_finish(iFrame)==0)
        MaskMB_display = currentImg*200;
        MaskMB_display(MaskMB_new>0)=255;
        
        figure(1);
        imagesc(MaskMB_display);  axis off; axis image; colormap(gray);
        hold on; plot(user_center_x(iFrame),user_center_y(iFrame),'r.');
        
        figure(2);
        hold off; imagesc(currentImg);  axis off; axis image; colormap(gray);
        RoiYX = bwboundaries(MaskMB_new);
        contourYX = RoiYX{1};
        hold on;
        plot(contourYX(:,2),contourYX(:,1),'r');
        % Wait for the user to input the cutting line, as the sequence of
        % points in the ROI input
        h = impoly;
        position = wait(h);
        if(size(position,1)<2)
            % if the user input just 1 points, means no more cutting is needed
            manual_finish(iFrame)=1;
        else
            % Cutting by defining all the pixels along this cutting line as
            % background
            for iP = 1 : size(position,1)-1                
               try
                   x1 = position(iP,1);
                y1 = position(iP,2);
                x2 = position(iP+1,1);
                y2 = position(iP+1,2);
                
                step_y = (y2-y1)/2/round((sqrt((y2-y1)^2+(x2-x1)^2)));
                step_x = (x2-x1)/2/round((sqrt((y2-y1)^2+(x2-x1)^2)));
                
               plot_x = x1:step_x:x2;
               plot_y = y1:step_y:y2;
               
               kill_connection_ind = sub2ind([height_image width_image],round(plot_y),round(plot_x));
               MaskMB_new(kill_connection_ind)=0;
               end
            end
            % Get rid of the "non-touching" parts, which is newly
            % disconnected regions due to wiping of the cutting line pixels
            labelMask = bwlabel(MaskMB_new,4);
            MaskMB_new = labelMask==labelMask(round(y),round(x));
            labelMask = bwlabel(MaskMB_new,4);
            obAreas = regionprops(labelMask,'Area');            
            if length(obAreas) > 1
                obAreas = [obAreas.Area];
                [dummy,iSort] = sort(obAreas,'descend');
                MaskMB_new = labelMask == iSort(1);
            end            
            MaskMB_new = imfill(MaskMB_new,'hole');                        
        end
    end
    
    % Save the resulting image to the hard disk
    imwrite(MaskMB_new,[Cell_manual_seg_folder,'/one_cells_open_',num2str(iFrame),'.TIF']);
    
end
  
                
  