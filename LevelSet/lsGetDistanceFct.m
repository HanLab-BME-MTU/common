function [val, val_matrix] = lsGetDistanceFct(mask_img, grid_coordinates, domain, x_spline, y_spline, signed)

% this function set the zero level curve, interpolation onto the 
% grid points based on the 2-d spline representing initial solution


if 0
   % Get the intersection of the x-gridlines 
   [x_x_intersection, x_y_intersection]= lsIntersectMaskLine(domain.x_grid_lines, mask_img, 1);

   % Improve the localization of the intersections
   [x_x_intersection, x_y_intersection]= lsIntersectSplineLine(domain.x_grid_lines, mask_img, 1);
   
   % Get the intersection of the y-gridlines 
   [y_x_intersection, y_y_intersection] = lsIntersectMaskLine(domain.y_grid_lines, mask_img, 2);


   % all intersection points
   x_intersection = cat(2, x_x_intersection, y_x_intersection);
   y_intersection = cat(2, x_y_intersection, y_y_intersection);

   figure,plot(x_intersection, y_intersection,'.');
   hold on
   p=1:x_spline.knots(end);
   x_spline_points = fnval(x_spline, p);
   y_spline_points = fnval(y_spline, p);
   plot(x_spline_points, y_spline_points,'r');
end


if 1
   % get the intersection of the x-gridlines with the edge
   x1(1:length(domain.x_grid_lines)) = 1;
   x2(1:length(domain.x_grid_lines)) = x_spline.knots(end);
   [x_approx_intersection_parameter, x_intersecting_lines] = lsIntersectApproxSplineLine(domain.x_grid_lines, x_spline, x1, x2);
   
   x1 = x_approx_intersection_parameter - 1;
   x2 = x_approx_intersection_parameter + 1;
   x_intersection_parameter = lsIntersectSplineLine(x_intersecting_lines, x_spline, x1, x2);

   % get the intersection of the y-gridlines with the edge
   y1(1:length(domain.y_grid_lines)) = 1;
   y2(1:length(domain.y_grid_lines)) = y_spline.knots(end);
   [y_approx_intersection_parameter, y_intersecting_lines] = lsIntersectApproxSplineLine(domain.y_grid_lines, y_spline, y1, y2);
   
   y1 = y_approx_intersection_parameter - 1;
   y2 = y_approx_intersection_parameter + 1;
   y_intersection_parameter = lsIntersectSplineLine(y_intersecting_lines, y_spline, y1, y2);

   x_x_intersection = fnval(x_spline,x_intersection_parameter);
   y_x_intersection = fnval(y_spline,x_intersection_parameter);

   x_y_intersection = fnval(x_spline,y_intersection_parameter);
   y_y_intersection = fnval(y_spline,y_intersection_parameter);

   % all intersection points
   x_intersection = cat(2, x_x_intersection, x_y_intersection);
   y_intersection = cat(2, y_x_intersection, y_y_intersection);

   figure,plot(x_x_intersection, y_x_intersection,'.');
   hold on
   plot(x_y_intersection,y_y_intersection,'.g');
   p=1:x_spline.knots(end);
   x_spline_points = fnval(x_spline, p);
   y_spline_points = fnval(y_spline, p);
   plot(x_spline_points, y_spline_points,'r');
end

% plot the grid lines
for i=1:length(domain.x_grid_lines)
   line([domain.x_grid_lines(i) domain.x_grid_lines(i)],[0 domain.y_size]);
end
for i=1:length(domain.y_grid_lines)
   line([0 domain.x_size], [domain.y_grid_lines(i) domain.y_grid_lines(i)],'Color','c');
end


% caliculate the minimal distance for each grid point
h_waitbar = waitbar(0,'Processing');
num_grid_coordinates = size(grid_coordinates,1);
for j = 1:num_grid_coordinates
    waitbar(j/num_grid_coordinates, h_waitbar, num2str(j));
    for i = 1:length(x_intersection)
        dist(i) = sqrt((grid_coordinates(j,1) - x_intersection(i))^2 + (grid_coordinates(j,2) - y_intersection(i))^2);
    end
    if mask_img(grid_coordinates(j,2), grid_coordinates(j,1)) > 0 & signed
       level_set_sign = -1; 
    else
       level_set_sign = 1;
    end
    val(j) = level_set_sign * min(dist);
end


% put the minimal distances from vetor into matrix form
val_matrix = reshape(val',size(domain.y_grid_lines,2),size(domain.x_grid_lines,2));




%surface(grid_coordinates(:,1),grid_coordinates(:,2),z);
%D   =  createDistanceMatrix(grid_coordinates(1,:), [x_intersection', y_intersection']);
 
% M and N are the matrices containing the set of point coordinates.
% M and N can represent point positions in 1, 2 and 3D, as follows.
%   
% 
% M=[ y1 x1     and   N=[ y1 x1
% y2 x2              y2 x2
% ...                ...
% ym xm ]            yn xn ]
% OUTPUT   D : distance matrix D=(dij), i=1..m, j=1..n