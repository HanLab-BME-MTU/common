function phi_zero = lsGetZeroLevel(phi, domain)
% LSGETZEROLEVEL finds the zero level countour in a 2D matrix
% 
%
% SYNOPSIS      phi_zero = lsGetZeroLevel(phi, domain)  
%
% INPUT         phi     : 2D matrix
%               domain  : structure with x and y grid line coordinates          
% 
% OUTPUT        phi_zero: zero level contour, 2D vector (x,y)    
%                     
%                           
% DEPENDENCES    lsLineMatching uses { contourc
%                                    } 
%
% Matthias Machacek 6/22/04

% this function gets the level at the grid lines ->
% x-levels, y-levels.
% it finds them by linear interpolation of the level set 
% matrix
% interpolate Level set

% interpolate on finer grid
%phi_fine = interp2(domain.x_grid_lines, domain.y_grid_lines', phi, domain.x_grid_lines_f, domain.y_grid_lines_f','cubic');

% extract zero level set
%phi_zero = contourc(domain.x_grid_lines_f, domain.y_grid_lines_f, phi_fine,[0 0]);
phi_zero = contourc(domain.x_grid_lines, domain.y_grid_lines, phi,[0 0]);

if ~isempty(phi_zero)
    phi_zero(:,1)=[];
end


% contourslice(...,[cv cv])

% alternative method:
if 0
   x_spacing = 0.1;
   y_spacing = 0.1;
   
   
   x_fine_grid_lines = 1:x_spacing:x_size;
   y_fine_grid_lines = 1:y_spacing:y_size;
   
   phi_fine = interp2(x_grid_lines, y_grid_lines, phi, x_fine_grid_lines, y_fine_grid_lines', 'spline');
   
   zero_level_approx = (phi_fine > -0.3) & (phi_fine < 0.3);
   figure
   imshow(zero_level_approx,[]);
   hold on
   contour(x_fine_grid_lines,y_fine_grid_lines, phi_fine, 40);
   
   
   figure,
   plot(zero_level_approx2(1,:),zero_level_approx2(2,:),'.'); 
   
   figure,
   plot(zero_level_approx2(1,:),zero_level_approx2(2,:),'.');
   
   % figure
   % contour(x_grid_lines,y_grid_lines, phi, 40);
end



