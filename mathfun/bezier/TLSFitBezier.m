function [P t r] = TLSFitBezier(X, n, varargin) %#ok<INUSD,STOUT>
% TLSFitBezier computes the best total least squares fit of a nth-degree
% Bezier curve to a set of ordered data points.
%
% Inputs:
% X       a m x d array representing a set of d-dimensional points
% n       the degree of the Bezier curve.
%
% Outputs:
% P       a n+1 x d array representing the set of d-dimensional control
%         points defining the Bezier curve.
%
% t       a m x 1 array of parameter value for each input points. t belongs
%         to [0, 1].
%
% res     a m x 1 array of residue, i.e the distance between each point and
%         the fitted curve.