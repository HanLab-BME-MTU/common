function track_points = lsIntegrateVelocity(dist_matrix, velocity_fct_matrix, grid_coordinates, delta_t, delta_x, delta_y, i_end, j_end, domain)


% find curve of interest (zero level)
phi_zero_t0 = lsGetZeroLevel(dist_matrix(:,:,1), domain); 

% discretize curve -> points of intereset
% take the points directly from lsGetZeroLevel function
track_points(:,:,1) = phi_zero_t0;

% find a B-spline interpolation of the velocity field
velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, velocity_fct_matrix(:,:,1));
fnplt(velocity_fct_spline), axis equal, axis off



for i = 1:size(dist_matrix,3)
    % find a B-spline interpolation of the velocity field
    velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, velocity_fct_matrix(:,:,i));
    
    % get speed at these points 
    track_points_velocity = fnval(velocity_fct_spline, track_points(:,:,i));
    
    figure 
    fnplt(velocity_fct_spline), axis equal
    %hold on
    %plot3(track_points(1,:,i), track_points(2,:,i), track_points_velocity, 'ro');
    
    % get the velocity direction at these points (grad phi)
    [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(dist_matrix(:,:,i), delta_x, delta_y, i_end, j_end);
    
    % find a B-spline interpolation of the gradient field
    grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x);
    grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y);
    
%     figure
%     fnplt(grad_x_spline)
%     hold on
%     fnplt(grad_y_spline)
    
    % get the gradient at the track points
    track_points_grad_x = fnval(grad_x_spline, track_points(:,:,i));
    track_points_grad_y = fnval(grad_y_spline, track_points(:,:,i));
    
   % figure
   % quiver(track_points(1,:,i), track_points(2,:,i), track_points_grad_x, track_points_grad_y);
    
    grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
    %track_points_grad_x = track_points_grad_x ./ grad;
    %track_points_grad_y = track_points_grad_y ./ grad;
    
    %delta_t*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j))
    
    % integrate velocity
    track_points(1,:,i+1) = track_points(1,:,i) - track_points_grad_x .* track_points_velocity * delta_t(i);
    track_points(2,:,i+1) = track_points(2,:,i) - track_points_grad_y .* track_points_velocity * delta_t(i);
end



