function [vertices,edges,edgePaths] = skel2graph(skelIn,nConn)
%SKEL2GRAPH converts a binary 3D skeleton matrix into a graph structure with nodes and edges 
% 
% [vertices,edges] = skel2graph(skelIn)
%
% [vertices,edges,edgePaths] = skel2graph(skelIn)
% 
% [vertices,edges,edgePaths] = skel2graph(skelIn,nConn)
% 
% Input:
% 
%   skelIn - The binary, 3D matrix containing a skeleton (maximum width is
%   1)
% 
%   nConn - Optional. The connectivity number to use (6,18 or 26)
% 
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
% 
% Hunter Elliott
% 4/2/2010


%% ---------- Input ----------- %%

%NOTE: At some point this should be generalized to 2D or 3D. Also, the
%exact paths the edges follow should be returned.

showPlots = true;

if nargin < 2 || isempty(nConn)
    nConn = 26;
end


%Validate that the input mask is in fact a skeleton(not too thick?)??


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

%Combine the end-points and junctions and label them
[vertMat,nVerts] = bwlabeln(vertMat | junctionPoints,nConn);

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

if nargin > 2
    edgePaths = cell(1,nEdges);
end

%Go through each edge...
for j = 1:nEdges
            
    %Find vertices which this edge connects
    edges(j,:) = unique(vertMat(imdilate(edgeMat == j,nHood) & vertMat));           
    
    %If requested, return the coordinates of each point on this edge
%     if nargout > 2
%         Get the coord of one vertex 
%         
%         
%         THIS DOESN"T SUPPORT 3D MATRICES - NEED TO WRITE OWN FUNC
%         
%         edgePaths{j} = bwtraceboundary(edgeMat == j & vertMat == edges(j,1),,'N');
%     end
end




if showPlots
   
    fsFigure(.5);
    title('Edges')    
    hold on
    spy3d(skelIn,'.k')
    cols = jet(nEdges);    
    arrayfun(@(x)(spy3d(edgeMat == x & skelIn,'.','color',cols(x,:),'MarkerSize',15)),1:nEdges)
    
    
    fsFigure(.5);
    title('Vertices')
    hold on
    spy3d(skelIn,'.k')
    cols = jet(nVerts);    
    arrayfun(@(x)(spy3d(vertMat == x & skelIn,'.','color',cols(x,:),'MarkerSize',15)),1:nVerts)
    
    
end

