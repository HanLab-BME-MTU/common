function x_velocity = ls_dis_interpolate(F, x, domain)
% interpolate F at the the positions x

num_points = size(x,1);


delta_x_inv = 1/domain.x_spacing;
delta_y_inv = 1/domain.y_spacing;

for i = 1:num_points
    % find the corners belonging to x(i)
    j=1;
    while x(i,1) >= domain.x_grid_lines(j)
        j=j+1;
    end
    
    x_2 = j;
    x_1 = j-1;
    
    j=1;
    while x(i,2) >= domain.y_grid_lines(j)
        j=j+1;
    end
    
    y_2 = j;
    y_1 = j-1;    
    
    
    % check for evolving edges in the x direction
    Dx_minus(x_1,y_1) = (F(x_1  ,y_1) - F(x_1-1y_1,)) * delta_x_inv;
    Dx_plus(x_2,y_1)  = (F(x_1+1,y_1) - F(x_1  ,y_1)) * delta_x_inv;
    
    
    % check for evolving edges in the y direction 
    Dy_minus(x_1,y_1) = (F(x_1,y_1  ) - F(x_1,y_1-1)) * delta_y_inv;
    Dy_plus(x_2,y_1)  = (F(x_1,y_1+1) - F(x_1,y_1  )) * delta_y_inv;    
    

    if Dx_minus(x_1,y_1) * Dx_plus(x_2,y_1) < 0
           % edge detected        
           % do a linear, one sided interpolation
           x_velocity(i) =   
        
        
    end
        
    
    
    
end
