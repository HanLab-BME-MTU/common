function [idx, dist] = KDTreeRangeQuery(inPts,queryPts,ranges) 
%KDTREERANGEQUERY finds all of the points which are within the specified radius of the query points
% 
% [idx, dist] = KDTreeRangeQuery(inPts,queryPts,ranges)
% 
% This function returns the indices of the input points which are within
% the specified range of the query points. Supports 2D or 3D point sets. In
% other words, this returns all the indices of the input points which are
% contained in the cuboids whose centers are given by the query points and
% whose dimensions are given by the input range vector.
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
%     ranges - an NxK matrix specifying the ranges from each query point to
%     find input points, e.g. the dimensions of the cuboids whithin which input
%     points will be found.
%     NOTE: This value should be of class double, or strange behavior may
%     occur.
% 
% 
% Output:
% 
%   idx - Nx1 cell array, the n-th element of which gives the indices of
%   the input points which are within the n-th range the n-th query
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
idx = cell(numQueryPoints, 1);
dist = cell(numQueryPoints, 1);
for i = 1:numQueryPoints
    curQueryRange = [queryPts(i,:) - 0.5 * ranges(i,:); queryPts(i,:) + 0.5 * ranges(i,:)];
    [idx{i}] = kdtreeobj.range( curQueryRange' );
    numNeighPoints = numel( idx{i} );
    ptNeigh = inPts( idx{i}, : );
    dist{i} = sqrt( sum( (ptNeigh - repmat(queryPts(i,:), numNeighPoints, 1)).^2, 2 ) );
end
