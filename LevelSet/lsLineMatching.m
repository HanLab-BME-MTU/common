function lsLineMatching
% LSLINEMATCHING matches two lines based on the Level Set method
% 
%
% SYNOPSIS       lsLineMatching
%
% INPUT                    
% 
% OUTPUT              
%                     
%                           
% DEPENDENCES    lsLineMatching uses { 
%                                    } 
%
% Matthias Machacek 6/10/04

test = 0;
if 1
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Test data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 200;
   domain.y_size = 200;
   
   domain.x_spacing = 15; 
   domain.y_spacing = 15;
   
   %create circle
   circle   = rsmak('circle',50,[0, 0]);
   ellipse  = fncmb(circle,[1.5 0;0 0.75]);
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
else
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   domain.x_size = 439;
   domain.y_size = 345;
   
   domain.x_spacing = 10; 
   domain.y_spacing = 10;
   
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
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the pixel edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   n_pix     = fscanf(fid_pixel_edge,'%g  ', [1 1]);
   x_p_edge  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix]);
   y_p_edge  = fscanf(fid_pixel_edge,'%g  ', [1 n_pix]);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   time = 1;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% read the mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   fileName_mask=char(filelist_mask(time));
   mask_img_t0=imread(fileName_mask);
   
   fileName_mask=char(filelist_mask(time+1));
   mask_img_t1=imread(fileName_mask);   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % End data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   x_spline_t0 = edge_sp_array_x(1);
   y_spline_t0 = edge_sp_array_y(1);
   
   x_spline_t1 = edge_sp_array_x(2);
   y_spline_t1 = edge_sp_array_y(2);
end



[grid_coordinates, x_grid, y_grid] = lsGenerateGrid(domain);

% fill the grid_line field in structure domain
domain = lsGenerateGridLines(domain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get edge at former time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signed_t0 = 1;
[val_t0, val_matrix_t0] = lsGetDistanceFct(mask_img_t0, grid_coordinates, domain, x_spline_t0, y_spline_t0, signed_t0);

% test the ddistance function and the zero level set finding
if test
    phi_zero = lsGetZeroLevelSet(val_matrix_t0, domain);
    figure,plot(phi_zero(1,:), phi_zero(2,:),'.');
    hold on
    plot(circle_points(1,:),circle_points(2,:),'r.'); 
    title('Original front, red, and extracted front, blue');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get edge at present time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signed_t1 = 1;
[val_t1, val_matrix_t1] = lsGetDistanceFct(mask_img_t1, grid_coordinates, domain, x_spline_t1, y_spline_t1, signed_t1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% plot results at previous time step
p_t0=1:x_spline_t0.knots(end);
x_spline_points_t0 = fnval(x_spline_t0, p_t0);
y_spline_points_t0 = fnval(y_spline_t0, p_t0);

figure,surface(domain.x_grid_lines,domain.y_grid_lines,val_matrix_t0);
hold on
plot(x_spline_points_t0, y_spline_points_t0,'r','LineWidth',2);
contour(domain.x_grid_lines,domain.y_grid_lines,val_matrix_t0, 40);


% plot results at present time step
p_t1=1:x_spline_t1.knots(end);
x_spline_points_t1 = fnval(x_spline_t1, p_t1);
y_spline_points_t1 = fnval(y_spline_t1, p_t1);

figure,surface(domain.x_grid_lines,domain.y_grid_lines,val_matrix_t1);
hold on
plot(x_spline_points_t1, y_spline_points_t1,'r','LineWidth',2);
contour(domain.x_grid_lines,domain.y_grid_lines,val_matrix_t1, 40);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Advance in time             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t = 0.05
delta_x = domain.x_spacing;
delta_y = domain.x_spacing;
i_end   = size(val_matrix_t1,1);
j_end   = size(val_matrix_t1,2);

phi_next = val_matrix_t0;

h_waitbar = waitbar(0,'Processing');
figure
hold on
plot(x_spline_points_t1, y_spline_points_t1,'r');
num_time_steps = 200;
for t=1:num_time_steps
   waitbar(t/num_time_steps, h_waitbar, num2str(t));
   [phi_next, delta_t_opt(t)] = lsSolveConvection(phi_next, delta_t, delta_x, delta_y, i_end, j_end, val_matrix_t1, domain);
   if mod(t,20) == 0 || t == 1
      phi_zero = lsGetZeroLevelSet(phi_next, domain);
      plot(phi_zero(1,:),phi_zero(2,:));
   end
end

figure
plot(delta_t_opt);
title('Optimal time step based on CFL number');
figure
contour(domain.x_grid_lines, domain.y_grid_lines, val_matrix_t0, 40);
hold on
contour(domain.x_grid_lines, domain.y_grid_lines, phi_next, 40);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 











