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
TEST_CASE_1 = 0;

if TEST_CASE_1
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
else
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
   time_increment = 1;
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




if 0
    % get the trail points
    trial_grid_points = lsFindTrailPoints(x_X_i_t0, y_X_i_t0, x_Y_i_t0, y_Y_i_t0, domain);

    % get the distance fct values of the trail points
    trial_grid_points(:,3) = lsGetDistanceFctVec(mask_img_t0,...
        trial_grid_points, known_zero_level_points_t0, domain, 1);

    % get dist_fct_tmp (Fast Marching Method)
    num_x_grid_lines = length(domain.x_grid_lines);
    num_y_grid_lines = length(domain.y_grid_lines);
    LARGE_NUMBER = 100000;
    dist_fct_tmp0(1:num_x_grid_lines,1:num_y_grid_lines) = LARGE_NUMBER;
    
    % get distance function by fast marching
    dist_fct_tmp = lsFastMarching(dist_fct_tmp0, trial_grid_points, domain, LARGE_NUMBER);
    figure
    surface(domain.x_grid_lines,domain.y_grid_lines,dist_fct_tmp);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get edge at former time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signed_t0 = 1;
[val_t0, val_matrix_t0] = lsGetDistanceFct(mask_img_t0, grid_coordinates,...
                            known_zero_level_points_t0, domain, signed_t0);

phi_zero_t0 = lsGetZeroLevel(val_matrix_t0, domain); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get distance function for all grid points at present time step %%
signed_t1 = 1;
[val_t1, val_matrix_t1] = lsGetDistanceFct(mask_img_t1, grid_coordinates,...
                             known_zero_level_points_t1, domain, signed_t1);

phi_zero_t1 = lsGetZeroLevel(val_matrix_t1, domain); 
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

phi(:,:,1) = val_matrix_t0;

h_waitbar = waitbar(0,'Processing');
figure
hold on
plot(x_spline_points_t1, y_spline_points_t1,'r');
max_time_steps = 100;
time_step = 1;
residual(time_step) = 100;
solution_difference(time_step) = 100;
c=jet(max_time_steps);

if 1
    % matlab time solver
    % Put into vector
    val_matrix_t0_vec = reshape(val_matrix_t0, prod(size(val_matrix_t0)),1);
    %options = odeset('OutputFcn',@odeplot, 'Events',@events);
    options = odeset('Events',@events);
    odeset;
    % 
    % ode45
    % stiff:
    % ode15s
    [t_steps, Y, TE,YE,IE] = ode45(@dy_fct,[0 20],val_matrix_t0_vec, options,...
        val_matrix_t1, i_end, j_end, delta_x, delta_y, domain);
    if isempty(TE)
        display('Solution did not converge in given time limit');
    else
        display(['Solution converged at time step  ' num2str(TE)]);
    end
    Y=Y';
    
    % Number of time steps
    num_time_steps = length(t_steps);
    for t=1:num_time_steps
        phi(:,:,t) = reshape(Y(:,t),i_end, j_end);

        % Re-calculate the velocities
        [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi(:,:,t),...
                delta_x, delta_y, i_end, j_end);
            
        velocity_fct(:,:,t) = lsGetVelocityFct(phi(:,:,t), val_matrix_t1,...
            i_end, j_end, grad_x, grad_y, domain);
    end



    % Extract the zero level
    phi_zero = lsGetZeroLevel(phi(:,:,end), domain);
    figure
    plot(phi_zero(1,:), phi_zero(2,:),'g');

    delta_t_opt = diff(t_steps);
    figure
    plot(delta_t_opt);
    title('Time steps');
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while time_step <= max_time_steps & solution_difference(time_step) > 0.1
        waitbar(time_step/max_time_steps, h_waitbar, num2str(time_step));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if time_step == 1
            [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(val_matrix_t0,...
                delta_x, delta_y, i_end, j_end);
            velocity_fct_0 = lsGetVelocityFct(val_matrix_t0, val_matrix_t1,...
                i_end, j_end, grad_x, grad_y, domain);

            [phi(:,:,time_step+1), velocity_fct(:,:,time_step), delta_plus, delta_minus, delta_t_opt(time_step)] = lsSolveConvection(phi(:,:,time_step),...
                phi(:,:,time_step), velocity_fct_0,...
                delta_plus, delta_minus,...
                delta_t, delta_x, delta_y, i_end, j_end, val_matrix_t1, domain);
        else
            [phi(:,:,time_step+1), velocity_fct(:,:,time_step), delta_plus, delta_minus, delta_t_opt(time_step)] = lsSolveConvection(phi(:,:,time_step),...
                phi(:,:,time_step-1), velocity_fct(:,:,time_step-1),...
                delta_plus, delta_minus,...
                delta_t, delta_x, delta_y, i_end, j_end, val_matrix_t1, domain);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extract the zero level
        phi_zero = lsGetZeroLevel(phi(:,:,time_step+1), domain);

        if mod(time_step,1) == 0 || time_step == 1
            plot(phi_zero(1,:),phi_zero(2,:),'Color',[c(time_step,1) c(time_step,2) c(time_step,3)]);
        end
        time_step = time_step+1;
        residual(time_step)            = norm(phi(:,:,time_step) - val_matrix_t1, 'fro');
        solution_difference(time_step) = norm(phi(:,:,time_step) - phi(:,:,time_step-1), 'fro');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    residual(1) = residual(2);
    plot(residual)
    title('Residual');

    figure
    solution_difference(1) = solution_difference(2);
    plot(solution_difference);
    title('Solution difference');
    
    figure
    plot(delta_t_opt);
    title('Optimal time step based on CFL number');  
end
time_step



figure
contour(domain.x_grid_lines, domain.y_grid_lines, val_matrix_t0, 40);
hold on
contour(domain.x_grid_lines, domain.y_grid_lines, phi(:,:,end), 40);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%% Integrate the velocity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
track_points = lsIntegrateVelocity(phi, velocity_fct, grid_coordinates,...
                      delta_t_opt, delta_x, delta_y, i_end, j_end, domain);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure
plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g');
hold on
plot(phi_zero(1,:), phi_zero(2,:),'r');
plot(phi_zero_t1(1,:), phi_zero_t1(2,:),'m');
for p = 1:size(track_points,2)
    plot(squeeze(track_points(1,p,:)), squeeze(track_points(2,p,:)), '-');   
end
axis equal
title('Tracked points');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure
plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g');
hold on
plot(phi_zero(1,:), phi_zero(2,:),'r');
for t = 1:size(track_points,3)
    plot(track_points(1,:,t), track_points(2,:,t), '-');   
end
if ~TEST_CASE_1
    plot(fnval(x_spline_tb,1: x_spline_tb.knots(end)),...
         fnval(y_spline_tb,1: y_spline_tb.knots(end)) ,'y');   
end
axis equal
title('Time lines');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy_vec = dy_fct(t, y, val_matrix_t1, i_end, j_end, delta_x, delta_y, domain)

phi_t = reshape(y, i_end, j_end);
[delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi_t, delta_x, delta_y, i_end, j_end);
F = lsGetVelocityFct(phi_t, val_matrix_t1, i_end, j_end, grad_x, grad_y, domain);

dy = zeros(i_end, j_end);
for i=1:i_end
    for j=1:j_end
        dy(i,j) = -max(F(i,j), 0) * delta_plus(i,j) -...
                   min(F(i,j), 0) * delta_minus(i,j);
    end
end
dy_vec = reshape(dy, numel(dy),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = events(t,y, phi_t1, i_end, j_end, delta_x, delta_y, domain)
% value(i) is the value of the function. 
% isterminal(i) = 1 if the integration is to terminate at a zero of this event 
%     function and 0 otherwise.
% direction(i) = 0 if all zeros are to be computed (the default), +1 if only 
%        the zeros where the event function increases, and -1 if only the 
%        zeros where the event function decreases.   

phi_t = reshape(y, i_end, j_end);
residual = norm(phi_t - phi_t1, 'fro');
if residual < 0.5
    value = 0;
else
    value = residual;
end
isterminal = 1;

direction = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = outputfcn(t,y, flag, phi_t1, i_end, j_end, delta_x, delta_y, domain)
status = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

