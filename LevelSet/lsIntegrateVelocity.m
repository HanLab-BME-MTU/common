function track_points = lsIntegrateVelocity(dist_matrix, velocity_fct_matrix, grid_coordinates, delta_t, delta_x, delta_y, i_end, j_end, domain)
% LSINTEGRATEVELOCITY integrates velocity to get time dependent position
%    
%
%
% SYNOPSIS   track_points = lsIntegrateVelocity(dist_matrix, velocity_fct_matrix, grid_coordinates, delta_t, delta_x, delta_y, i_end, j_end, domain)
%
%
% INPUT      dist_matrix            :
%            velocity_fct_matrix    :
%            grid_coordinates       :
%            delta_t                :   
%            delta_x                :
%            delta_y                :
%            i_end                  :
%            j_end                  :
%            domain                 :
%                          
% 
% OUTPUT     track_points           :
%              
%                           
% DEPENDENCES    lsIntegrateVelocity uses {                                
%                                       }
%
%                lsIntegrateVelocity is used by { 
%                                           }
%
% Matthias Machacek 06/24/04


track_points(:,:,1) = lsGetZeroLevel(dist_matrix(:,:,1), domain);


for i = 1:size(dist_matrix,3)-1
    % find a B-spline interpolation of the velocity field
    velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, velocity_fct_matrix(:,:,i));
    
    % get speed at these points 
    track_points_velocity = fnval(velocity_fct_spline, track_points(:,:,i));
    
    %figure 
    %track_points_velocityfnplt(velocity_fct_spline), axis equal
    %hold on
    %plot3(track_points(1,:,i), track_points(2,:,i), track_points_velocity, 'ro');
    av_vel(i) = mean(track_points_velocity);
    
    % get the velocity direction at these points (grad phi)
    [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(dist_matrix(:,:,i), delta_x, delta_y, i_end, j_end);
    
    % find a B-spline interpolation of the gradient field
    grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x);
    grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y);
    
%     figure
%     quiver(domain.x_grid_lines, domain.y_grid_lines, grad_x, grad_y);
%     hold on
%     plot(track_points(1,:,i), track_points(2,:,i)); 
%     axis equal
%     quiver(domain.y_grid_lines, domain.x_grid_lines, grad_y, grad_x);
    
    delta_plus_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, delta_plus);
    delta_minus_spline = csapi({domain.x_grid_lines, domain.y_grid_lines},delta_minus);
    track_points_delta_plus =   fnval(delta_plus_spline, track_points(:,:,i));
    track_points_delta_minus =  fnval(delta_minus_spline, track_points(:,:,i));
    
%     figure
%     fnplt(grad_x_spline)
%     hold on
%     fnplt(grad_y_spline)
    
    % get the gradient at the track points
    track_points_grad_x = fnval(grad_x_spline, track_points(:,:,i));
    track_points_grad_y = fnval(grad_y_spline, track_points(:,:,i));
    
    %figure
    %quiver(track_points(1,:,i), track_points(2,:,i), track_points_grad_x, track_points_grad_y);
    
    grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
    track_points_grad_x_u = track_points_grad_x ./ grad;
    track_points_grad_y_u = track_points_grad_y ./ grad;
    
    %delta_t*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j))
    
    % integrate velocity
    track_points(1,:,i+1) = track_points(1,:,i) -  track_points_grad_x_u .* track_points_velocity * delta_t(i);
    track_points(2,:,i+1) = track_points(2,:,i) -  track_points_grad_y_u .* track_points_velocity * delta_t(i);
end



