function dist = distancePointBezier(cP,P,t0,t1)
% distancePointBezier computes the distance between a point and a linear,
% quadratic or cubic 3D Bezier curve. By specifying a parametrization 
% interval with the optional inputs the distance to a segment
% can be computed.
%
% Required Inputs:
% cP                A N x 3 array representing a set of 3-dimensional 
%                   control points where N is 2,3 or 4          
% 
% Optional Inputs:
% t0 = 0 (default)  Start value of the parametrization interval            
% t1 = 1 (default)  End value of the parametrization interval
%
% Outputs:
% dist              The distance to the Bezier curve
%
% Pascal Berard, May 2011