function [s_approx, lines_red] = lsIntersectApproxSplineLine(lines, spline, x1, x2)
% LSINTERSECTSPLINELINE finds all approximately intersections between set of lines and spline 
%    
%
%
% SYNOPSIS      intersection_parameter = lsIntersectApproxSplineLine(lines, spline, x1, x2)
%
% INPUT         lines  :  array with line coordinate
%               spline :  one component spline
%               x1     :  lower spline parameter boundary
%               x2     :  upper spline parameter boundary
% 
% OUTPUT        intersection_parameter : parameter of the pline at the
%                                        intersecions
%                           
% DEPENDENCES   lsIntersectSplineLine uses {   fminbnd                             
%                                       }
%
%               lsIntersectSplineLine is used by { lsGetDistanceFct
%                                           }
%
% Matthias Machacek 06/09/04



% find first the approximate intersection


s_approx = [];
lines_red = [];
for i = 1:length(lines)
   s = x1(i):0.2:x2(i);
   f = (lines(i) - fnval(spline,s)).^2;
   [min_dist, min_index] = find(f < 1); 
   s_approx = cat(2,s_approx,s(min_index));
   if length(min_index) > 0
      for j=1:length(min_index)
         lines_red(end+1) = lines(i);
      end
   end
end