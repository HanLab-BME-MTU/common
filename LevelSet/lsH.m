function dy_vec = lsH(t, y, val_matrix_t1, i_end, j_end, delta_x, delta_y, domain)

if 1
    phi_vec = y(1:i_end*j_end);
    phi_t = reshape(phi_vec , i_end, j_end);
    x_vec = y(i_end*j_end+1 : end);
    x = reshape(x_vec, length(x_vec)/2,2);
else
    phi_t = reshape(y, i_end, j_end);
end

[delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi_t, delta_x, delta_y, i_end, j_end);
%[delta_plus, delta_minus, grad_x, grad_y] = weno_5(phi_t, delta_x, delta_y, i_end, j_end);
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
    %velocity_fct_spline = spapi({6 6}, {domain.x_grid_lines, domain.y_grid_lines}, F');
   
    % Get velocity at these points
    x_velocity = fnval(velocity_fct_spline, x');    
    %x_velocity = interp2(domain.x_grid_lines,domain.y_grid_lines, F ,x(:,1), x(:,2),'cubic');
    %x_velocity = x_velocity';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use the matlab build-in gradient function
    %[grad_x, grad_y] = gradient(phi_t, delta_x, delta_y);
    [grad_x, grad_y] = lsGradient(phi_t, 4, 0, delta_x, delta_y, i_end, j_end);

    % Find a B-spline interpolation of the gradient field
    grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x');
    grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y');
    % grad_x_spline = spapi({6 6}, {domain.x_grid_lines, domain.y_grid_lines}, grad_x');
    % grad_y_spline = spapi({6 6}, {domain.x_grid_lines, domain.y_grid_lines}, grad_y');

    % Get the gradient at the track points
     track_points_grad_x = fnval(grad_x_spline, x');
     track_points_grad_y = fnval(grad_y_spline, x');
%    track_points_grad_x = interp2(domain.x_grid_lines,domain.y_grid_lines, grad_x ,x(:,1), x(:,2),'cubic');
%    track_points_grad_y = interp2(domain.x_grid_lines,domain.y_grid_lines, grad_y ,x(:,1), x(:,2),'cubic');
%    track_points_grad_x = track_points_grad_x';
%    track_points_grad_y = track_points_grad_y';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
    
    dx_x = x_velocity .* track_points_grad_x./ grad;
    dx_y = x_velocity .* track_points_grad_y./ grad;
    
    dy_vec = cat(1, dy_vec, dx_x');
    dy_vec = cat(1, dy_vec, dx_y');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%