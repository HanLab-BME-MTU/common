function len = lengthBezier(P,t0,t1)
% function len = lengthBezier(P,t0,t1)
% lengthBezier computes the length of a 3D Bezier curve. By specifying a 
% parametrization interval with the optional inputs the length of a segment
% can be computed.
%
% Required Inputs:
% P                 A N x 3 array representing a set of 3-dimensional control points
% 
% Optional Inputs:
% t0 = 0 (default)  Start value of the parametrization interval            
% t1 = 1 (default)  End value of the parametrization interval
%
% Outputs:
% len               The length of the Bezier curve
%
% Pascal Berard, May 2011