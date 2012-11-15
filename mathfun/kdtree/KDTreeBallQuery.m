function [idx, dist] = KDTreeBallQuery(inPts,queryPts,radii)
%KDTREEBALLQUERY finds all of the points which are within the specified radius of the query points
% 
% [idx, dist] = KDTreeBallQuery(inPts,queryPts,radii)
% 
% This function returns the indices of the input points which are within
% the specified radii of the query points. Supports 1D, 2D or 3D point sets. 
% In other words, this returns all the indices of the input points which are
% contained in the spheres whose centers are given by the query points and
% whose radii are given by the input radii vector.
%
% Input:
% 
%     inPts - an MxK matrix specifying the input points to test for distance
%     from the query points, where M is the number of points and K is the
%     dimensionality of the points.
% 
%     queryPts - an NxK matrix specifying the query points, e.g. the centers of
%     the spheres within which input points will be found.
% 
%     radii - an Nx1 vector or a scalar specifying the distances from each 
%     query point to find input points, e.g. the radii of the spheres within 
%     which input points will be found. If scalar, the same radius is used
%     for all query points.
%     NOTE: This value should be of class double, or strange behavior may
%     occur.
% 
% 
% Output:
% 
%   idx - Nx1 cell array, the n-th element of which gives the indices of
%   the input points which are within the n-th radii of the n-th query
%   point.
% 
%   dist - Nx1 cell array, the n-th element of which gives the corresponding 
%   distances between the input points and the n-th query point.
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
%        O(n^(1-1/d) + m) time where m is the number of reported points
%        and d is the dimensionality of the points.   
% 
%     So you would be acting contrary to the purpose (fast-querying) of using KDTree 
%     if you are going to call this function multiple times for the same KD-tree. 
%     Use the KDTree class in KDTree.m instead.
%       

warning( ['The KD-Tree is built each time you call this function.' ...
         'If you want to build the kdtree once and query it multiple times, then ' ...
         'it would be really suboptimal to use this function because building a KD-tree ' ...
         'takes a lot more time than querying it. Use the KDTree class in KDTree.m instead. ' ...
         'We may be removing this function from the toolkit in the near future.'] );
         
kdtreeobj = KDTree( inPts );

numQueryPoints = size(queryPts,1);
if isscalar( radii )
    radii = radii + zeros(numQueryPoints,1);
end

idx = cell(numQueryPoints, 1);
dist = cell(numQueryPoints, 1);
for i = 1:numQueryPoints
    [idx{i}, dist{i}] = kdtreeobj.ball( queryPts(i,:), radii(i) );
end

