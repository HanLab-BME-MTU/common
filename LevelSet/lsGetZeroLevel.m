function phi_zero = lsGetZeroLevelSet(phi, domain)

phi_zero = contourc(domain.x_grid_lines, domain.y_grid_lines, phi,[0 0]);
phi_zero(:,1)=[];

% figure,
% plot(phi_zero(1,:), phi_zero(2,:), '.');
% axis equal





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



