function [cP,t] = truncateBezier(cP,tStart,tEnd,t)
% function [cP,t] = truncateBezier(cP,tStart,tEnd,t)
% truncateBezier computes the control points of a segment of the input 
% Bezier curve. Optionally, points of the old curve can be mapped onto the
% new curve. In comparison to truncateBezierOld, the input curve can be of 
% any degree, however it might be a bit slower. The implementation is based on
% the SISL library.
%
% Required Inputs:
% cP                A N x D array representing a set of 3-dimensional 
%                   control points where N is 2,3 or 4. D is the dimension
%                   of the control points.
% tStart            Defines the start point of the segment            
% tEnd              Defines the end point of the segment
% 
% Optional Inputs:
% t                 M x 1 array containing parameter values of the old
%                   curve that will be mapped to the new one
%
% Outputs:
% cP                The control points of the new curve
% t                 The transformed parametrization 
%
% Pascal Berard, March 2012

