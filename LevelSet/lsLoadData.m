function [mask_img_t0, mask_img_t1,x_spline_t0, y_spline_t0, x_spline_t1,y_spline_t1,...  
    known_zero_level_points_t0, known_zero_level_points_t1, grid_coordinates, domain] = lsLoadData(TEST_CASE)


if TEST_CASE == 1
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = 5; 
   domain.y_spacing = 5;
   
   %create circle
   circle   = rsmak('circle',50,[0, 0]);
   ellipse  = fncmb(circle,[1.5 0;0 0.6]);
   circle   = fncmb(circle,[0.6 0;0 1.5]);
   %ellipse  = fncmb(circle,[1.5 0;0 1.5]);
   circle   = fncmb(circle,'+', 100);
   ellipse  = fncmb(ellipse,'+',100);
   fnplt(circle);
   hold on
   fnplt(ellipse);
   axis equal
   
   % Create mask
   p = 0:0.05:circle.pieces;
   circle_points = fnval(circle,p);
   p = 0:0.05:ellipse.pieces;
   ellipse_points = fnval(ellipse,p);
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, circle_points(1,:)', circle_points(2,:)');
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, ellipse_points(1,:)', ellipse_points(2,:)');

   % create x,y splines
   s_p = 1:length(p);
   x_spline_t0 = fn2fm(spline(s_p, circle_points(1,:)),'B-');
   y_spline_t0 = fn2fm(spline(s_p, circle_points(2,:)),'B-');
   
   x_spline_t1 = fn2fm(spline(s_p, ellipse_points(1,:)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, ellipse_points(2,:)),'B-');
   
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0] = lsGetGridIntersections(x_spline_t0, y_spline_t0, domain);
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif TEST_CASE == 2
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = 5; 
   domain.y_spacing = 5;
   
   %create circle
   circle   = rsmak('circle',62,[0, 0]);
   ellipse  = fncmb(circle,[1.5 0;0 0.6]);
   ellipse  = fncmb(ellipse,'+', [100]);
   
   circle1   = rsmak('circle',30,[0, 0]);
   circle2   = rsmak('circle',30,[0, 0]);
   circle1   = fncmb(circle1,'+', [65 100]);
   circle2   = fncmb(circle2,'+', [135 100]);
   
   fnplt(circle1);
   hold on
   fnplt(circle2);   
   fnplt(ellipse);
   axis equal
   
   % Create mask
   p = 0:0.05:circle1.pieces;
   circle_points1 = fnval(circle1,p);
   p = 0:0.05:circle2.pieces;
   circle_points2 = fnval(circle2,p);   
   circle_points = cat(2, circle_points1,circle_points2); 
   p = 0:0.05:ellipse.pieces;
   ellipse_points = fnval(ellipse,p);
   
   
   mask_img_t0 = roipoly(domain.y_size, domain.x_size, circle_points(1,:)', circle_points(2,:)');
   mask_img_t1 = roipoly(domain.y_size, domain.x_size, ellipse_points(1,:)', ellipse_points(2,:)');

   % create x,y splines
   s_p = 1:length(p);
   x_spline_t01 = fn2fm(spline(s_p, circle_points1(1,:)),'B-');
   y_spline_t01 = fn2fm(spline(s_p, circle_points1(2,:)),'B-');
   
   x_spline_t02 = fn2fm(spline(s_p, circle_points2(1,:)),'B-');
   y_spline_t02 = fn2fm(spline(s_p, circle_points2(2,:)),'B-');   
   
   x_spline_t1 = fn2fm(spline(s_p, ellipse_points(1,:)),'B-');
   y_spline_t1 = fn2fm(spline(s_p, ellipse_points(2,:)),'B-');
   
   x_spline_t0 = fncmb(x_spline_t01,'+',x_spline_t02);
   y_spline_t0 = fncmb(y_spline_t01,'+',y_spline_t02);
   x_spline_t0 = fn2fm(x_spline_t0,'B-');
   y_spline_t0 = fn2fm(y_spline_t0,'B-');
   % End test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   % (find the known points)
   [x_X_i_t01, y_X_i_t01, x_Y_i_t01, y_Y_i_t01] = lsGetGridIntersections(x_spline_t01, y_spline_t01, domain);
   [x_X_i_t02, y_X_i_t02, x_Y_i_t02, y_Y_i_t02] = lsGetGridIntersections(x_spline_t02, y_spline_t02, domain);
   x_X_i_t0 = cat(2,x_X_i_t01, x_X_i_t02);
   y_X_i_t0 = cat(2,y_X_i_t01, y_X_i_t02);
   x_Y_i_t0 = cat(2,x_Y_i_t01, x_Y_i_t02);
   y_Y_i_t0 = cat(2,y_Y_i_t01, y_Y_i_t02);
   
   known_zero_level_points_t0(:,1) = [x_X_i_t0'; x_Y_i_t0'];
   known_zero_level_points_t0(:,2) = [y_X_i_t0'; y_Y_i_t0'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get intersection of grid with the curve %%%%%%%%%%%%%%%%%%%%%%%%%
   [x_X_i_t1, y_X_i_t1, x_Y_i_t1, y_Y_i_t1] = lsGetGridIntersections(x_spline_t1, y_spline_t1, domain);
   known_zero_level_points_t1(:,1) = [x_X_i_t1'; x_Y_i_t1'];
   known_zero_level_points_t1(:,2) = [y_X_i_t1'; y_Y_i_t1'];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif TEST_CASE == 3
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 439;
   domain.y_size = 345;
   
   domain.x_spacing = 5; 
   domain.y_spacing = 5;
   
   %cd L:\projects\rho_protrusion\cell1\protrusion
   cd /lccb/projects/rho_protrusion/cell1/protrusion_1-30_ps2_s40
   load edge_spline
   
   PROJECT_DIR = '/lccb/projects/rho_protrusion/cell1/';
   PROT_DIR = 'protrusion_1-30_ps2_s40/';
   IMG_NAME = 'plane01';
   
   % file name containing the edge pixels
   file_pixel_edge=[PROJECT_DIR  PROT_DIR 'pixel_edge.dat'];
   
   fid_pixel_edge= fopen(file_pixel_edge,'r');
   if fid_pixel_edge == -1
      error('Could not open file pixel edge');
   end
   
   % mask file names
   firstfilename_mask=([PROJECT_DIR PROT_DIR 'mask_' IMG_NAME '.tif']);
   [filelist_mask]=getFileStackNames(firstfilename_mask);
   
   time = 4;
   time_increment = 10;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1 : time
       n_pix_t0     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       y_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       known_zero_level_points_t0 = [x_p_edge_t0', y_p_edge_t0'];
   end
   for i = 1 : time_increment
       n_pix_t1     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       y_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       known_zero_level_points_t1 = [x_p_edge_t1', y_p_edge_t1'];
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+time_increment));
   mask_img_t1=imread(fileName_mask);   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(time);
   y_spline_t0 = edge_sp_array_y(time);
  
   x_spline_tb = edge_sp_array_x(time+round(time_increment/2));
   y_spline_tb = edge_sp_array_y(time+round(time_increment/2));
   
   x_spline_t1 = edge_sp_array_x(time+time_increment);
   y_spline_t1 = edge_sp_array_y(time+time_increment);
   
   
    [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
    
    
    
elseif TEST_CASE == 4
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 439;
   domain.y_size = 345;
   
   domain.x_spacing = 5; 
   domain.y_spacing = 5;
   
   cd /lccb/projects/alpha/W512r_smaller/protrusion_01-90_s30_p20
   load edge_spline
   
   PROJECT_DIR = '/lccb/projects/alpha/W512r_smaller/';
   PROT_DIR = 'protrusion_01-90_s30_p20/';
   IMG_NAME = 'cut_cut_cut_W512r01';
   
   % file name containing the edge pixels
   file_pixel_edge=[PROJECT_DIR  PROT_DIR 'pixel_edge.dat'];
   
   fid_pixel_edge= fopen(file_pixel_edge,'r');
   if fid_pixel_edge == -1
      error('Could not open file pixel edge');
   end
   
   % mask file names
   firstfilename_mask=([PROJECT_DIR PROT_DIR 'mask_' IMG_NAME '.tif']);
   [filelist_mask]=getFileStackNames(firstfilename_mask);
   
   time = 4;
   time_increment = 10;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1 : time
       n_pix_t0     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       y_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       known_zero_level_points_t0 = [x_p_edge_t0', y_p_edge_t0'];
   end
   for i = 1 : time_increment
       n_pix_t1     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       y_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       known_zero_level_points_t1 = [x_p_edge_t1', y_p_edge_t1'];
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+time_increment));
   mask_img_t1=imread(fileName_mask);   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(time);
   y_spline_t0 = edge_sp_array_y(time);
  
   x_spline_tb = edge_sp_array_x(time+round(time_increment/2));
   y_spline_tb = edge_sp_array_y(time+round(time_increment/2));
   
   x_spline_t1 = edge_sp_array_x(time+time_increment);
   y_spline_t1 = edge_sp_array_y(time+time_increment);
   
   
    [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);
    
elseif TEST_CASE == 5
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

   
   cd /lccb/projects/alpha/PtK1_control/cut_s399/protrusion_1-243_s40_p20
   load edge_spline
   
   PROJECT_DIR = '/lccb/projects/alpha/PtK1_control/cut_s399/';
   PROT_DIR = 'protrusion_1-243_s40_p20/';
   IMG_NAME = 'cut_S399ACTIN001';
   
   % Read image to determine the size
   c_image = imread('img_edge_cut_S399ACTIN001.tif');
   [img_h, img_w, c] = size(c_image);
   
   domain.x_size = img_w;
   domain.y_size = img_h;
   
   domain.x_spacing = 5; 
   domain.y_spacing = 5;
   
   
   % file name containing the edge pixels
   file_pixel_edge=[PROJECT_DIR  PROT_DIR 'pixel_edge.dat'];
   
   fid_pixel_edge= fopen(file_pixel_edge,'r');
   if fid_pixel_edge == -1
      error('Could not open file pixel edge');
   end
   
   % mask file names
   firstfilename_mask=([PROJECT_DIR PROT_DIR 'mask_' IMG_NAME '.tif']);
   [filelist_mask]=getFileStackNames(firstfilename_mask);
   
   time = 4;
   time_increment = 10;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for i = 1 : time
       n_pix_t0     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       y_p_edge_t0  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t0]);
       known_zero_level_points_t0 = [x_p_edge_t0', y_p_edge_t0'];
   end
   for i = 1 : time_increment
       n_pix_t1     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
       x_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       y_p_edge_t1  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix_t1]);
       known_zero_level_points_t1 = [x_p_edge_t1', y_p_edge_t1'];
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+time_increment));
   mask_img_t1=imread(fileName_mask);   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(time);
   y_spline_t0 = edge_sp_array_y(time);
  
   x_spline_tb = edge_sp_array_x(time+round(time_increment/2));
   y_spline_tb = edge_sp_array_y(time+round(time_increment/2));
   
   x_spline_t1 = edge_sp_array_x(time+time_increment);
   y_spline_t1 = edge_sp_array_y(time+time_increment);
   
   
    [grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

    % fill the grid_line field in structure domain
    domain = lsGenerateGridLines(domain);  
end

