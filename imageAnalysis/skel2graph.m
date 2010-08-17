function [vertices,edges] = skel2graph(skelIn,nConn)
%SKEL2GRAPH converts a binary 3D skeleton matrix into a graph structure with nodes and edges 
% 
% [vertices,edges] = skel2graph(skelIn)
%
%                        ... = skel2graph(skelIn,nConn)
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


%   edgePaths - An Nx1 cell array containing the ordered coordinates of
%   each point along each edge. NOTE: This will make the processing MUCH
%   slower!
% 
% Hunter Elliott
% 4/2/2010


%% ---------- Input ----------- %%

%NOTE: At some point this should be generalized to 2D or 3D. Also, the
%exact paths the edges follow should be returned.

showPlots = false;

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

%Label end-points first. This is done separately because it allows
%recognition of "spurs" as endpoints - branches which are 1 voxel in
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

% if nargout > 2
%     edgePaths = cell(nEdges,1);
%     edgeInit = zeros(3*max(size(skelIn)),3); %Matrix for over-initializing edge paths    
% end

%Go through each edge...
for j = 1:nEdges
            
    %Find vertices which this edge connects
    tmp = unique(vertMat(imdilate(edgeMat == j,nHood) & vertMat));           
    if length(tmp) == 2
        edges(j,:) = tmp;
    elseif length(tmp) > 2
        error('Problem with input matrix! Check that it is in fact a skeleton, and that it''s connectivity  matches the specified connectivity!')
    end
      
%COMMENTED THIS OUT BECAUSE ITS REALLY SLOW AND FAILS IF THE SKELETON HAS LOOPS    
%     %If requested, return the coordinates of each point on this edge. I'm
%     %pretty sure there's a faster way to do this...???
%     if nargout > 2
%         %Get the coord of one vertex        
%         iVert = 1;
%         edgePaths{j} = edgeInit;        
%         currPos = vertices(edges(j,1),:);
%         while 1            
%             tmpMask = false(size(skelIn));
%             tmpMask(round(currPos(1)),round(currPos(2)),round(currPos(3))) = true;
%             iNeighbors = find(imdilate(tmpMask,nHood) & (edgeMat == j));
%             if length(iNeighbors) == 1
%                 edgePaths{j}(iVert,:) = currPos;
%                 iVert = iVert + 1;
%                 [currPos(1),currPos(2),currPos(3)] = ind2sub(size(skelIn),iNeighbors);                
%             else
%                 break
%             end
%         end
%         edgePaths{j} = edgePaths{j}(1:iVert-1,:);
%    end

end




if showPlots
   
    fsFigure(.5);
    title('Edges are points, vertices are circles')    
    hold on    
    cols = lines(nEdges);    
    arrayfun(@(x)(spy3d(edgeMat == x & skelIn,'.','color',cols(x,:),'MarkerSize',5)),1:nEdges)            
    cols = lines(nVerts);
    arrayfun(@(x)(spy3d(vertMat == x & skelIn,'o','color',cols(x,:),'MarkerSize',5)),1:nVerts)    
    
end

