function [dist,t] = distancePointBezierSISL(cP,P)
% function [dist,t] = distancePointBezierSISL(cP,P)
% distancePointBezierSISL computes the distance between a point and a Bezier
% curve. If only one control point is specified the distance between the 
% point and the control point is returned. In comparison to distancePointBezier
% the parametrization interval cannot be specified. However, there is no limitation
% on the complexity of the curves. This function uses the SISL NURBS
% library which can be found at "http://www.sintef.no/sisl".
%
% Required Inputs:
% cP                A N x 3 array representing a set of 3-dimensional 
%                   control points where N is 1,2,3 or 4          
% 
% Optional Inputs:
%
% Outputs:
% dist              The distance to the Bezier curve
% t                 The corresponding curve parameter
%
% Pascal Berard, March 2012