function [phi_t1, phi_t1_calc, track_points, protrusion] = lsLineMatching(mask_img_t0, mask_img_t1, x_spline_t0, y_spline_t0, x_spline_t1, y_spline_t1)
% LSLINEMATCHING matches two lines based on the Level Set method
% 
%
% SYNOPSIS       lsLineMatching
%
% INPUT                    
% 
% OUTPUT         phi_t1         : given distance function at second time step
%                phi(:,:,end)   : calculated distance function at second time step
%                track_points   : points on the edge (x,y,t)
%
% DEPENDENCES    lsLineMatching uses { 
%                                    } 
%
% Matthias Machacek 6/10/04

CONTROL = 1; 

%TEST_CASE = 1; % two ellipses
%TEST_CASE = 2; % two non-intersecting circles
%TEST_CASE = 3; % ellipses two cirlces
 TEST_CASE = 4; % whole cell 
%TEST_CASE = 5; % part cell 
%TEST_CASE = 6; % part cell cut_s399
 
INT_VELOCITY = 1;
RESULT_DIR = '/lccb/projects/alpha/Level_set_test/';

[mask_img_t0, mask_img_t1, x_spline_t0, y_spline_t0, x_spline_t1,y_spline_t1,...  
    known_zero_level_points_t0, known_zero_level_points_t1, grid_coordinates, domain] = lsLoadData(TEST_CASE, CONTROL);


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
[val_t0, phi_t0] = lsGetDistanceFct(mask_img_t0, grid_coordinates,...
                            known_zero_level_points_t0, domain, signed_t0);

phi_zero_t0 = lsGetZeroLevel(phi_t0, domain); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get distance function for all grid points at present time step %%
signed_t1 = 1;
[val_t1, phi_t1] = lsGetDistanceFct(mask_img_t1, grid_coordinates,...
                             known_zero_level_points_t1, domain, signed_t1);

phi_zero_t1 = lsGetZeroLevel(phi_t1, domain); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%if TEST_CASE ~= 2
%     p_t01=1:x_spline_t01.knots(end);
%     x_spline_points_t01 = fnval(x_spline_t01, p_t01);
%     y_spline_points_t01 = fnval(y_spline_t01, p_t01);
%     p_t02=1:x_spline_t02.knots(end);
%     x_spline_points_t02 = fnval(x_spline_t02, p_t01);
%     y_spline_points_t02 = fnval(y_spline_t02, p_t01);
% 
%     x_spline_points_t0 = cat(2, x_spline_points_t01, x_spline_points_t02);
%     y_spline_points_t0 = cat(2, y_spline_points_t01, y_spline_points_t02);
% end

p_t0=1:x_spline_t0.knots(end);
x_spline_points_t0 = fnval(x_spline_t0, p_t0);
y_spline_points_t0 = fnval(y_spline_t0, p_t0);

p_t1=1:x_spline_t1.knots(end);
x_spline_points_t1 = fnval(x_spline_t1, p_t1);
y_spline_points_t1 = fnval(y_spline_t1, p_t1);

if CONTROL & 0
    % plot results at previous time step
    h_phi_t0 = figure;
    surface(domain.x_grid_lines,domain.y_grid_lines, phi_t0);
    hold on
    plot(x_spline_points_t0, y_spline_points_t0,'r','LineWidth',2);
    contour(domain.x_grid_lines,domain.y_grid_lines, phi_t0, 40);

    % plot results at present time step
    h_phi_t1 = figure;
    surface(domain.x_grid_lines,domain.y_grid_lines,phi_t1);
    hold on
    plot(x_spline_points_t1, y_spline_points_t1,'r','LineWidth',2);
    contour(domain.x_grid_lines,domain.y_grid_lines,phi_t1, 40);
    
    % plot the two contours
    h_edges = figure;
    plot(x_spline_points_t0, y_spline_points_t0, 'b');
    hold on
    plot(x_spline_points_t1, y_spline_points_t1, 'r');
        
    % Save  results
    hgsave(h_phi_t0, [RESULT_DIR 'phi_t0.fig']); 
    print(h_phi_t0,  [RESULT_DIR 'phi_t0.eps'],'-depsc2','-tiff'); 
    print(h_phi_t0,  [RESULT_DIR 'phi_t0.tif'],'-dtiff');
    hgsave(h_phi_t1, [RESULT_DIR 'phi_t0.fig']); 
    print(h_phi_t1,  [RESULT_DIR 'phi_t0.eps'],'-depsc2','-tiff'); 
    print(h_phi_t1,  [RESULT_DIR 'phi_t0.tif'],'-dtiff');  
    hgsave(h_edges,  [RESULT_DIR 'edges_t0_t1.fig']); 
    print(h_edges,   [RESULT_DIR 'edges_t0_t1.eps'],'-depsc2','-tiff'); 
    print(h_edges,   [RESULT_DIR 'edges_t0_t1.tif'],'-dtiff');    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Advance in time             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_x = domain.x_spacing;
delta_y = domain.x_spacing;
i_end   = size(phi_t1,1);
j_end   = size(phi_t1,2);

figure
plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'r');
hold on
   
% Put into vector
phi_t0_vec = reshape(phi_t0, prod(size(phi_t0)),1);  
if INT_VELOCITY
    track_points_0 = phi_zero_t0';
    
    %p_t0=1:0.4:x_spline_t0.knots(end);
    %track_points_0(:,1) = fnval(x_spline_t0, p_t0);
    %track_points_0(:,2) = fnval(y_spline_t0, p_t0);
    
    num_track_points = size(track_points_0, 1);
    phi_t0_vec = cat(1, phi_t0_vec,track_points_0(:,1));
    phi_t0_vec = cat(1, phi_t0_vec,track_points_0(:,2));
end

%options = odeset('RelTol',1e-4,'AbsTol',1e-7,'OutputFcn',@outputfcn,'Events',@events);
options = odeset('OutputFcn',@outputfcn, 'Events',@events);
%options = odeset('Events',@events);%'MaxStep',0.1,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t_steps, Y, TE,YE,IE] = ode45(@dy_fct,[0 200],phi_t0_vec, options,...
    phi_t1, i_end, j_end, delta_x, delta_y, domain);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(TE)
    display('Solution did not converge in given time limit');
else
    display(['Solution converged at time step  ' num2str(TE)]);
end
Y=Y';

% Number of time steps
num_time_steps = length(t_steps);

% Initialize the displacement
protrusion = zeros((size(Y,1) - (i_end * j_end))/2,  1);

for t=1:num_time_steps
    if INT_VELOCITY
        phi_vec = Y(1:i_end * j_end,t);
        phi(:,:,t) = reshape(phi_vec,i_end, j_end);
        track_points_vec = Y(i_end * j_end+1 : end,t);
        track_points(:,:,t) = reshape(track_points_vec, num_track_points, 2);
    else
        phi(:,:,t) = reshape(Y(:,t),i_end, j_end);
    end
end

% Initialize the displacement
protrusion = zeros((size(Y,1) - (i_end * j_end))/2,  1);

for t=1:num_time_steps - 1
    % get the displacements
    dx = track_points(:,1,t) - track_points(:,1,t+1);
    dy = track_points(:,2,t) - track_points(:,2,t+1);
    protrusion = protrusion + sqrt(dx .^ 2 + dy .^ 2); 
end
% Get the sign of the protrusion
%mask_img_t0f = flipud(mask_img_t0);
for i=1:size(track_points,1)
    if mask_img_t0(round(track_points(i,1,end)), round(track_points(i,2,end))) == 0
        prot_sign = 1;
    else
        prot_sign = -1;
    end 
    protrusion(i) = prot_sign .* protrusion(i);
end

phi_t1_calc = phi(:,:,end);
phi_zero = lsGetZeroLevel(phi_t1_calc, domain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CONTROL
    h_f = figure(gcf);
    hgsave(h_f, [RESULT_DIR 'zero_level_evol.fig']);
    print(h_f,  [RESULT_DIR 'zero_level_evol.eps'],'-depsc2','-tiff');
    print(h_f,  [RESULT_DIR 'zero_level_evol.tif'],'-dtiff');
    
     % Show the time steps
    delta_t_opt = diff(t_steps);
    figure
    plot(delta_t_opt);
    title('Time steps');
    
    % Show the curvature of the last distance function
    kappa = lsCurvature(phi(:,:,end), delta_x, delta_y, i_end, j_end);
    figure
    surface(kappa);
    
    % Superimposed contour plots
    h_contour = figure;
    contour(domain.x_grid_lines, domain.y_grid_lines, phi_t0, 40);
    hold on
    contour(domain.x_grid_lines, domain.y_grid_lines, phi(:,:,end), 40);
    plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g','LineWidth',2);
    plot(phi_zero_t1(1,:), phi_zero_t1(2,:),'r','LineWidth',2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_tracked_points = figure;
    plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g');
    hold on
    plot(phi_zero(1,:), phi_zero(2,:),'r');
    plot(phi_zero_t1(1,:), phi_zero_t1(2,:),'--m');
    plot(x_spline_points_t0, y_spline_points_t0,'r+', 'MarkerSize',3);
    plot(x_spline_points_t1, y_spline_points_t1,'r+', 'MarkerSize',3);
    for p = 1:size(track_points,1)
        plot(squeeze(track_points(p,1,:)), squeeze(track_points(p,2,:)), '-');
    end
    axis equal
    title('Tracked points');
    legend('Zero level t0', 'Zero level t1 solution', 'Zero level t1 given',...
        'org. curve t0','org. curve t1')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_time_lines = figure
    plot(phi_zero_t0(1,:), phi_zero_t0(2,:),'g');
    hold on
    plot(phi_zero(1,:), phi_zero(2,:),'r');
    for t = 1:size(track_points,3)
        plot(track_points(:,1,t), track_points(:,2,t), '-');
    end
    % if TEST_CASE == 3
    %     plot(fnval(x_spline_tb,1: x_spline_tb.knots(end)),...
    %          fnval(y_spline_tb,1: y_spline_tb.knots(end)) ,'y');
    % end
    axis equal
    title('Time lines');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_disp = figure;
    plot(protrusion);
    title('Protrusion [pixel/frame rate]');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save  results
    hgsave(h_tracked_points, [RESULT_DIR 'tracked_points.fig']);
    print(h_tracked_points,  [RESULT_DIR 'tracked_points.eps'],'-depsc2','-tiff');
    print(h_tracked_points,  [RESULT_DIR 'tracked_points.tif'],'-dtiff');
    hgsave(h_time_lines, [RESULT_DIR 'time_lines.fig']);
    print(h_time_lines,  [RESULT_DIR 'time_lines.eps'],'-depsc2','-tiff');
    print(h_time_lines,  [RESULT_DIR 'time_lines.tif'],'-dtiff');
    hgsave(h_disp, [RESULT_DIR 'protrusion.fig']);
    print(h_disp,  [RESULT_DIR 'protrusion.eps'],'-depsc2','-tiff');
    print(h_disp,  [RESULT_DIR 'protrusion.tif'],'-dtiff');  
    hgsave(h_contour, [RESULT_DIR 'contour.fig']);
    print(h_contour,  [RESULT_DIR 'contour.eps'],'-depsc2','-tiff');
    print(h_contour,  [RESULT_DIR 'contour.tif'],'-dtiff');      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy_vec = dy_fct(t, y, val_matrix_t1, i_end, j_end, delta_x, delta_y, domain)

if 1
    phi_vec = y(1:i_end*j_end);
    phi_t = reshape(phi_vec , i_end, j_end);
    x_vec = y(i_end*j_end+1 : end);
    x = reshape(x_vec, length(x_vec)/2,2);
else
    phi_t = reshape(y, i_end, j_end);
end

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

if 1
    % Integrate the velocity 
    %%%% Velocity interpolation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find a B-spline interpolation of the velocity field
    velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, F');

    % Get velocity at these points
    x_velocity = fnval(velocity_fct_spline, x');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use the matlab build-in gradient function
    %[grad_x, grad_y] = gradient(phi_t, delta_x, delta_y);
    [grad_x, grad_y] = lsGradient(phi_t, 10, 0, delta_x, delta_y, i_end, j_end);

    % Find a B-spline interpolation of the gradient field
    grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x');
    grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y');

    % Get the gradient at the track points
    track_points_grad_x = fnval(grad_x_spline, x');
    track_points_grad_y = fnval(grad_y_spline, x');

    grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dx_x = x_velocity .* track_points_grad_x./ grad;
    dx_y = x_velocity .* track_points_grad_y./ grad;
    
    dy_vec = cat(1, dy_vec, dx_x');
    dy_vec = cat(1, dy_vec, dx_y');
end
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

if 1
    phi_vec = y(1:i_end*j_end);
    phi_t = reshape(phi_vec , i_end, j_end);
    x_vec = y(i_end*j_end+1 : end);
    x = reshape(x_vec, length(x_vec)/2,2);
else
    phi_t = reshape(y, i_end, j_end);
end

residual = norm(phi_t - phi_t1, 'fro');
if residual < 1
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
f = figure(gcf);

if ~isempty(y)
    if 1
        phi_vec = y(1:i_end*j_end, end);
        phi = reshape(phi_vec , i_end, j_end);
        x_vec = y(i_end*j_end+1 : end, end);
        x = reshape(x_vec, length(x_vec)/2,2);
        % Extract the zero level
        phi_zero = lsGetZeroLevel(phi, domain);
        plot(phi_zero(1,:), phi_zero(2,:),'g');       
    else
        phi = reshape(y, i_end, j_end);
        % Extract the zero level
        phi_zero = lsGetZeroLevel(phi(:,:,end), domain);
        plot(phi_zero(1,:), phi_zero(2,:),'g');       
    end
end

status = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

