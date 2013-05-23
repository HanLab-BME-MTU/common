function [idx, dist] = KDTreeClosestPoint(inPts,queryPts)
%KDTREECLOSESTPOINT for every query point in queryPts, find the closest point belonging to inPts
% 
% [idx, dist] = KDTreeClosestPoint(inPts,queryPts)
% 
% This function returns the index of the input point closest to each inPts.
% Supports ND point sets.
%
% Input:
% 
%     inPts - an MxK matrix specifying the input points to test for distance
%     from the query points, where M is the number of points and K is the
%     dimensionality of the points.
% 
%     queryPts - an NxK matrix specifying the query points.
% 
% Output:
% 
%   idx - Nx1 array, the n-th element of which gives the index of
%   the input point closest to the the n-th query point.
% 
%   dist - Nx1 array, the n-th element of which gives the corresponding 
%   distance between the closest input point and the n-th query point.
%
% Warning:
% 
%     The KD-Tree is built each time you call this function. 
%     If you want to build the kdtree once and query it multiple times, then
%     it would be suboptimal to use this function because building a KD-tree 
%     takes a lot more time than querying it. Specifically,
%     
%     -- building the KD-Tree takes O(n log(n)) time
%     -- just range querying an already built tree takes on the average 
%        O( log(n) ).
% 
%     So you would be acting contrary to the purpose (fast-querying) of using KDTree 
%     if you are going to call this function multiple times for the same KD-tree. 
%     Use the KDTree class in KDTree.m instead.
%       

warning('KDTREE:closestPtBuild',['The KD-Tree is built each time you call this function.' ...
         'If you want to build the kdtree once and query it multiple times, then ' ...
         'it would be really suboptimal to use this function because building a KD-tree ' ...
         'takes a lot more time than querying it. Use the KDTree class in KDTree.m instead. ' ...
         'We may be removing this function from the toolkit in the near future.'] );
         
kdtreeobj = KDTree( inPts );
[idx, dist] = kdtreeobj.nn( queryPts );
