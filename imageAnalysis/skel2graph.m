function [vertices,edges,edgePaths] = skel2graph(skelIn,nConn)
%SKEL2GRAPH converts a binary 3D skeleton matrix into a graph structure with nodes and edges 
% 
% [vertices,edges] = skel2graph(skelIn)
% [vertices,edges,edgePaths] = skel2graph(skelIn)
%                        ... = skel2graph(skelIn,nConn)
% 
% Input:
% 
%   skelIn - The binary, 3D matrix containing a skeleton (maximum width is
%   1)
% 
%   nConn - Optional. The connectivity number to use (6,18 or 26). 
%           WARNING: I have only tested this function with 26-connected
%           skeletons (the kind produced by skeleton3D)!!
%           Optional. Default is 26.
% 
% Output:
% 
%   vertices - A Mx3 matrix of the coordinates of vertices in the
%   skeleton/graph (end-points or junction points), where M is the number
%   of vertices
% 
%   edges - An Nx2 matrix of the index of the vertices that each edge
%   connects.
%
%   edgePaths - An Nx1 cell array containing the ordered coordinates of
%   each point along each edge. NOTE: Requesting this output will make the
%   processing quite a bit slower!
% 
% Hunter Elliott
% 4/2/2010
%

%% ---------- Input ----------- %%

%NOTE: At some point this should be generalized to handle both 2D and
%3D.-HLE

showPlots = false;

if nargin < 1 || isempty(skelIn) || ~islogical(skelIn) || ndims(skelIn) ~= 3
    error('The first input must be a 3D binary matrix containing a skeleton!')
end

if nargin < 2 || isempty(nConn)
    nConn = 26;
end


%TEMP - should Validate that the input mask is in fact a skeleton(not too thick?)??-HLE


%% ------- Detection -----  %%
%Detect edges and vertices

%Get the neighborhood for this connectivity
nHood = bwnHood3D(nConn);

%Calculate local neighbor number
nNeighbors = bwNneighbors(skelIn,nHood);

%Find end-points
vertMat = (nNeighbors == 1) & skelIn;

%Find junction points
junctionPoints = (nNeighbors > 2) & skelIn;

%Label end-points first. This is done separately because it allows
%recognition of "spurs" as endpoints - "branches" which are 1 voxel in
%length. Otherwise these are combined with the adjacent vertex.
[tmp,nEP] = bwlabeln(vertMat,nConn);

%Label the junction points
[vertMat,nJP] = bwlabeln(junctionPoints,nConn);

%Combine them with the endpoints
vertMat(vertMat(:)>0) = vertMat(vertMat(:)>0) + nEP;%Shift the labels
vertMat(tmp(:)>0) = tmp(tmp(:)>0);%Add the endpoint labels. These override any junction-labels.
nVerts = nJP+nEP;

%Get edges
edgeMat = nNeighbors == 2 & skelIn;

%Label these edges
[edgeMat,nEdges] = bwlabeln(edgeMat,nConn);


%% ------ Connectivity ----- %%
%Determine which edges are connected to which vertices

edges = zeros(nEdges,2);
vertices = zeros(nVerts,3);

%Go through each vertex (these may in fact be clusters of points depending
%on the connectivity and skeleton structure)
for j = 1:nVerts
    
    %Get the index of the point(s)
    currInd = find(vertMat == j);
    
    %Convert this to matrix coord
    [currM,currN,currP] = ind2sub(size(skelIn),currInd);
    
    %Average/store these coord
    vertices(j,:) = [mean(currM),mean(currN),mean(currP)];        


end

if nargout > 2
    edgePaths = cell(nEdges,1);
    edgeInit = zeros(3*max(size(skelIn)),3); %Matrix for over-initializing edge paths    
end

%Go through each edge...
for j = 1:nEdges
            
    %Find vertices which this edge connects
    tmp = unique(vertMat(imdilate(edgeMat == j,nHood) & vertMat));           
    if length(tmp) == 2
        edges(j,:) = tmp;
    elseif length(tmp) > 2
        error('Problem with input matrix! Check that it is in fact a skeleton, and that it''s connectivity  matches the specified connectivity!')
    end
       
    %If requested, return the coordinates of each point on this edge. I'm
    %pretty sure there's a faster way to do this...???
    if nargout > 2
        %Get the coord of one vertex        
        iVert = 1;        
        
        %First, make sure it's not a spur
        if ~any(edges(j,:)==0)
            
            edgePaths{j} = edgeInit;

            %Try to find a vertex that is a single voxel and use it as a start
            %point for the edge path
            if nnz(vertMat == edges(j,1)) == 1
                currPos = round(vertices(edges(j,1),:));
            elseif nnz(vertMat == edges(j,2)) == 1
                currPos = round(vertices(edges(j,2),:));
            else
                %If both vertices are multi-voxel, we need to find an
                %appropriate start-point.

                %First, get coord of all edge points
                iEdgePts = find(edgeMat == j);
                edgeCoord = zeros(numel(iEdgePts),3);
                [edgeCoord(:,1),edgeCoord(:,2),edgeCoord(:,3)] = ind2sub(size(skelIn),iEdgePts);
                %Find the one closest to the first vertex
                [~,iClosest] = min(arrayfun(@(x)(sqrt(sum((edgeCoord(x,:) - ...
                                   vertices(edges(j,1),:)).^2))),1:numel(iEdgePts)));
                currPos = edgeCoord(iClosest,:);
            end
            %"Walk" along the edge from this start point until we hit the end
            iNeighbors = 1;
            while numel(iNeighbors) == 1            
                edgePaths{j}(iVert,:) = currPos;
                tmpMask = false(size(skelIn));
                tmpMask(round(currPos(1)),round(currPos(2)),round(currPos(3))) = true;
                iNeighbors = find(imdilate(tmpMask,nHood) & (edgeMat == j));            
                if numel(iNeighbors) == 1
                    iVert = iVert + 1;
                    %Delete the previous point to avoid loops/direction changes
                    edgeMat(currPos(1),currPos(2),currPos(3)) = 0;                                
                    [currPos(1),currPos(2),currPos(3)] = ind2sub(size(skelIn),iNeighbors);                
                end
            end
            %Remove extra points from over-initialization.
            edgePaths{j} = edgePaths{j}(1:iVert,:);
        end
    end
   
end



if showPlots %shows plots of outputs, for debugging
   
    fsFigure(.5); %#ok<UNRCH>
    cols = jet(nEdges);    
    hold on    
    if nargout > 2
        title('Edges are lines, vertices are circles')    
        arrayfun(@(x)(plot3(edgePaths{x}(:,2),edgePaths{x}(:,1),edgePaths{x}(:,3),'color',cols(x,:))),1:nEdges)            
    else
        title('Edges are points, vertices are circles')    
        arrayfun(@(x)(spy3d(edgeMat == x & skelIn,'.','color',cols(x,:),'MarkerSize',5)),1:nEdges)            
    end
    cols = lines(nVerts);
    arrayfun(@(x)(spy3d(vertMat == x & skelIn,'o','color',cols(x,:),'MarkerSize',5)),1:nVerts)    
    
end

