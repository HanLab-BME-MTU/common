function c_point = getClosestPt(point, sp_x, sp_y)

SOLUTION_SPACING = 0.5;

% generate points
r_n_lower = sp_x.knots(1);
r_n_upper = sp_x.knots(end);      
r_n  = r_n_lower: SOLUTION_SPACING : r_n_upper;


curve_points = [fnval(sp_x,r_n); fnval(sp_y,r_n)]';


% get distances
dist_cand = sqrt((curve_points(:,1) - point(1)).^2+ (curve_points(:,1) - point(2)).^2);
 
%get the minimal distance
[dist min_index] = min(dist_cand);

%get the point
c_point = curve_points(min_index,:);
