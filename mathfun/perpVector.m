function v_perp=perpVector(point1,e_1,point2)
%finds the perpendicular vector v connecting the straight line defined by
%point1 and e_1 with point2
%
%SYNOPSIS v_perp=perpVector(point1,e_1,point2)
%
%INPUT point1, e_1 : point and UNIT vector defining the straight line
%      point2: point from which the new perpendicular vector should start
%      (all inputs can be lists)
%
%OUTPUT v_perp perpendicular vector from the straight line to the point
%           (norm(v_perp)=distance of the point from the line)
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%e_1 has to be a list of unit vectors!

%vector between point1 and point2
v_temp=point2-point1;

%projection of v_temp on e_1
v_proj=(sum(v_temp.*e_1,2)*ones(1,3)).*e_1;

%v_temp - the parallel projection = the perpendicular part
v_perp=v_temp-v_proj;
