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

contr = 0;

track_points(:,:,1) = lsGetZeroLevel(dist_matrix(:,:,1), domain);

% number of time points 
num_time_steps = size(dist_matrix,3)-1;

for i = 1:num_time_steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Velocity interpolation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find a B-spline interpolation of the velocity field
    velocity_fct_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, velocity_fct_matrix(:,:,i)');
    
    % Get velocity at these points 
    track_points_velocity = fnval(velocity_fct_spline, track_points(:,:,i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Gradient field interpolation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the gradient at these points (grad phi)
    [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(dist_matrix(:,:,i), delta_x, delta_y, i_end, j_end);    
    
    % Find a B-spline interpolation of the gradient field
    grad_x_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_x');
    grad_y_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, grad_y');
    
    % Get the gradient at the track points
    track_points_grad_x = fnval(grad_x_spline, track_points(:,:,i));
    
    
    track_points_grad_y = fnval(grad_y_spline, track_points(:,:,i));
    	if contr
        figure 
        fnplt(velocity_fct_spline), axis equal
        hold on
        plot3(track_points(1,:,i), track_points(2,:,i), track_points_velocity, 'ro');
    end
    if contr
        [m_grad_x, m_grad_y] = gradient(dist_matrix(:,:,i), delta_x, delta_y);
        
        figure
        quiver(domain.x_grid_lines, domain.y_grid_lines, grad_x, grad_y,0);
        hold on
        contour(domain.x_grid_lines, domain.y_grid_lines, dist_matrix(:,:,i));
        %quiver(domain.x_grid_lines, domain.y_grid_lines, m_grad_x, m_grad_y);
        plot(track_points(1,:,i), track_points(2,:,i),'r');
        %axis equal
        %quiver(domain.y_grid_lines, domain.x_grid_lines, grad_y, grad_x);
        quiver(track_points(1,:,i), track_points(2,:,i), track_points_grad_x, track_points_grad_y,0,'r');
        axis equal
        xlabel('x');
        ylabel('y');
        
        for ii=1:size(grad_x,1)
            for jj=1:size(grad_x,2)
                grad_x_vec((ii-1)*size(grad_x,2)+jj) = grad_x(ii,jj);
                grad_y_vec((ii-1)*size(grad_y,2)+jj) = grad_y(ii,jj);
                x_cord((ii-1)*size(grad_x,2)+jj) = domain.x_grid_lines(ii);
                y_cord((ii-1)*size(grad_x,2)+jj) = domain.y_grid_lines(jj); 
            end
        end
         
        figure   
        fnplt(grad_x_spline)
        hold on
        plot3(track_points(1,:), track_points(2,:), track_points_grad_x,'o');
        plot3(x_cord, y_cord, grad_x_vec, 'x');
        xlabel('x');
        ylabel('y');
        
        figure
        fnplt(grad_y_spline)
        hold on
        plot3(track_points(1,:), track_points(2,:), track_points_grad_y,'o');
        plot3(x_cord, y_cord, grad_y_vec, 'x');
        xlabel('x');
        ylabel('y');
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     delta_plus_spline = csapi({domain.x_grid_lines, domain.y_grid_lines}, delta_plus);
%     delta_minus_spline = csapi({domain.x_grid_lines, domain.y_grid_lines},delta_minus);
%     track_points_delta_plus =   fnval(delta_plus_spline, track_points(:,:,i));
%     track_points_delta_minus =  fnval(delta_minus_spline, track_points(:,:,i));
    
    
    grad = sqrt(track_points_grad_x.^2 + track_points_grad_y.^2);
    track_points_grad_x_u = track_points_grad_x ./ grad;
    track_points_grad_y_u = track_points_grad_y ./ grad;
    
    %delta_t*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Integrate velocity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    track_points(1,:,i+1) = track_points(1,:,i) +  track_points_grad_x_u .* track_points_velocity * delta_t(i);
    track_points(2,:,i+1) = track_points(2,:,i) +  track_points_grad_y_u .* track_points_velocity * delta_t(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



