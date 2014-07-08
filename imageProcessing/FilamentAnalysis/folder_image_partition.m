% an example wrapper of the partition code

% Liya 06.30.2014

inImg = imread('peppers.png');
inImg = rgb2gray(inImg);

row_partition_number=4;
col_partition_number=5;

outImg_Cell = single_image_partition(inImg, row_partition_number, col_partition_number);

output_dir = 'C:/test';

figure(1);
imagesc(inImg); axis image; axis off; colormap(gray);

figure(2);

for iR = 1 : row_partition_number
    for iC = 1 : col_partition_number
        
        % the plotting is for debug testing
        subplot(row_partition_number,col_partition_number,sub2ind([col_partition_number,row_partition_number],iC,iR));
        imshow(outImg_Cell{iR, iC}); axis image; axis off; colormap(gray);
        
        %save the partitioned image to disk
        imwrite(outImg_Cell{iR, iC}, ...
            [output_dir,filesep,'part_R',num2str(iR),'_C',num2str(iC),'.tif']);
    end
end


% command name 
% copyfile
% movefile
% mkdir
% exist(output_dir, 'dir')

% h = hist(outImg_Cell{iR, iC}(:), 0:1:10);


