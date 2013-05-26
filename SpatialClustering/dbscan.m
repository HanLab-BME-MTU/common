%[clusterID] = dbscan(X, minsize, R) implements the DBSCAN clustering algorithm
%
% Inputs: 
%            X : Nx2 matrix of point coordinates
%      minsize : minimum accepted cluster size
%            R : neighborhood search radius
%
% Outputs:
%    clusterID : cluster identifier/class number for each point in X
%           nn : #neighbors within R for each point
%
%
% Adapted from the pseudo-code published in
% [1] M. Ester et al., "A density-based algorithm for discovering clusters in
% large spatial databases with noise," KDD-96 Proc., 1996.

% Francois Aguet, 05/26/2013


function [clusterID, nn] = dbscan(X, minsize, R)

np = size(X,1);
kdidx = KDTreeBallQuery(X, X, R);
nn = cellfun(@numel, kdidx)-1;

clusterID = NaN(np,1);

visited = false(np,1);

c = 1; % class/cluster label
for i = 1:np
    
    if ~visited(i)
        % neighbors of current point
        neighborPts = kdidx{i};
        
        % if point has < minsize neighbors -> outlier
        if numel(neighborPts) < minsize
            clusterID(i) = 0;
            visited(i) = true;
        else
            clusterID(neighborPts) = c;
            
            % loop through neighbors and add their neighbors to current cluster
            while ~isempty(neighborPts)
                % mark current neighbor as visited
                visited(neighborPts(1)) = true;
                
                % neighbors of current neighbor
                ni = kdidx{neighborPts(1)};
                                
                % remove current neighbor from list
                neighborPts(1) = [];
                
                if numel(ni) > 1
                    clusterID(ni) = c;
                    % add neighbors of current neighbor if not previously visited
                    nvisit = ni(~visited(ni));
                    neighborPts = [neighborPts; nvisit]; %#ok<AGROW>
                    visited(nvisit) = true;
                end
            end
            c = c+1;
        end
    end
end
