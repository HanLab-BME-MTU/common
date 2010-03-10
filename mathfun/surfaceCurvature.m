function K = surfaceCurvature(S,N)
%SURFACECURVATURE calculates the local curvature of each face in the input triangular mesh 
% 
% K = surfaceCurvature(surface,normals)
% 
% This function will calculate an approximate curvature value for each face
% in the input triangular mesh. The surface should follow the format used
% by patch, where the surface is contained in a structure with the fields
% "faces" and "veritces"
% 
% Input: 
% 
%   surface - The surface to calculate curvature on, using the FV
%   (Faces/vertices format) used by patch, isosurface etc.
% 
%   normals - The normals of the surface at each vertex. 
% 
%   Example:
%   To calculate the local curvature of an isosurface of an image, use the
%   following commands:
% 
%       s = isosurface(image,isoValue);
%
%       n = isonormals(image,s.vertices); 
% 
%       c = surfaceCurvature(s,n);
% 
%   Which can then be visualized with the command:
% 
%       patch(s,'FaceColor','flat','EdgeColor','none','FaceVertexCData',c)    
% 
% 
% Output:
% 
%   K = A Mx1, where M is the number of faces, vector of the approximate
%   gaussian curvature at each face.
% 
% 
% 
%Hunter Elliott 
%3/2010
%

if nargin < 2 || isempty(S) || isempty(N)
    error('Must input surface mesh and surface normals!')
end

%Number of faces
nTri = size(S.faces,1);

%Barycentric coordinates for location of interpolated normal
abc = ones(1,3) * 1/3; %This will estimate curvature at the center of each face.

%Init array for curvature values
K = zeros(nTri,1);

%Should probably vectorize this at some point...
for i = 1:nTri
    
    %Get the coordinates of this triangle's vertices
    X = S.vertices(S.faces(i,:),:);
    
    %Get the normal vectors for these vertices
    n = N(S.faces(i,:),:);
    
    %Interpolated normal
    ni = abc * n;
    
    %Triangle normal
    m = cross(X(2,:)-X(1,:),X(3,:)-X(2,:));
    
    %Gaussian curvature
    K(i) = det(n) / (dot(ni,ni)*dot(ni,m));
              
end


