function varargout = solvePDEVectorBoundary(xy,uv,pdePar,imgSize,meshQuality)
%SOLVEPDEVECTORBOUNDARY solves the selected PDE using the input vectors as a boundary condition 
% 
%                     [X,Y] = solvePDEVectorBoundary(xy,uv)
% [ux,px,ex,tx,uy,py,ey,ty] = solvePDEVectorBoundary(xy,uv)
%                     ...   = solvePDEVectorBoundary(xy,uv,pdePar)
%                     ...   = solvePDEVectorBoundary(xy,uv,pdePar,imSize)
%                     ...   = solvePDEVectorBoundary(xy,uv,pdePar,imSize,meshQuality)
% 
% This function solves the selected partial differential equation (PDE)
% over the area enclosed by the vectors x and y given the vector-valued
% boundary conditions specified by the vectors u and v. The X and Y
% components of the solution are solved separately as separate boundary
% conditions. The solution are found using the matlab PDE toolbox. The PDE
% must be of the form supported by the PDE toolbox function adaptmesh.m.
% 
% 
% Input:
% 
%   xy - The 2xM or Mx2 matrix containing the x and y positions of the
%         boundary of the area to solve the PDE. These positions must form
%         a closed polygon (M>=3 and the last and first points are
%         adjacent).
% 
%   uv - The 2xM or Mx2 vectors containing the x and y components of the
%         boundary condition at each boundary point specified by xy.
%         Must be the same length as xy.
% 
%   pdePar - 1x3 Vector specifying the parameters of the PDE to solve in
%            the object interior. The elements of this vector correspond to
%            the coefficients c, a and f used by the functions assembpde.m
%            and adaptmesh.m. Optional. Default is [1 0 0];
%
%   imgSize - A 1x2 positive integer vector containing the size of the area
%             over which to solve the PDE. Areas outside of the area
%             specified by x and y will have zero values. Only used if X,Y
%             outputs are requested. Optional. If not specified, the X and
%             Y matrices will be just large enough to fit area enclosed by
%             xy.
%
%   meshQualtiy - Integer scalar between 1-10. The quality of the
%                 triangular mesh to use when solving the PDE. Large
%                 numbers will DRASTICALLY increase computation time, while
%                 decreasing error in the solution.
%
%
% Output:
%
%   [X,Y] - The rectangular matrices specifying the X and Y components of
%           the solution to the PDE at regularly spaced points within the
%           solution area. This output is slightly slower as it requires
%           that the solution defined on the triangular mesh be
%           interpolated to a homogeneous grid. 
%
%      ------------------------ OR -----------------------------
%
%   [ux,px,ex,tx,uy,py,ey,ty] - The solution to the PDE on the triangular
%           mesh used by the PDE toolbox. This contains the solution values
%           u, vertex locations p, edge topology e and triangles t for the
%           x and y components of the solution.
%
% Hunter Elliott
% 8/2010
%

%% ----------- Input ---------- %%
if nargin < 2 || isempty(xy) || isempty(uv)
    error('Must input boundary location xy and boundary condition uv!')
end

%Convert to Mx2 if input as 2xM
if size(xy,2) ~= 2
    xy = xy';
end
if size(uv,2) ~= 2
    uv = uv';
end

%Check that border is closed. 
if ~isequal(xy(1,:),xy(end,:))
    error('The input polygon defined by xy must be closed - that is, the first and last points must be identical!')
end

nPts = size(xy,1);

if nPts < 3 || size(uv,1) ~= nPts
    %NOTE: It is not actually required that xy and uv be the same length,
    %but it is an easy way to ensure correspondence between a location on
    %the boundary and the boundary condition.... HLE
    error('xy and uv must be matrices of the same size - Mx2 or 2xM where M >= 3!')
end

if nargin < 3 || isempty(pdePar)
    pdePar = [1 0 0 ];
elseif numel(pdePar) ~= 3
    error('The input PDE parameters must be a 1x3 vector!')
end
    
if nargin < 4
    imgSize = [];
end

if nargin < 5 || isempty(meshQuality)
    meshQuality = 3;
elseif (meshQuality <= 0) || (round(abs(meshQuality)) ~= meshQuality) ...
        || meshQuality > 10
    error('Mesh quality must be a positive integer between 1 and 10!')
end

%% --------- Init ----------- %%    

%We use global variables for the boundary coord and values, because the PDE
%toolbox is very strict about geometry and boundary condition inputs.
global BOUND_COND
global OBJ_BOUND

%Check handedness of boundary and reverse if necessary. The boundary
%function assumes correct handedness when determining object outside vs
%inside


%Fit interpolating spline to boundary.
OBJ_BOUND = spline(linspace(0,1,nPts),xy');

%If the image size wasn't input, make sure the range includes the whole
%cell
if isempty(imgSize)
    sTmp = linspace(0,1,nPts);
    vs = ppval(OBJ_BOUND,sTmp);
    imgSize = ceil(max(vs,[],2) + 10); %Leave a little room...
end

%% -------- Solve ---------- %%
%Solves X and Y components of eqtn seperately


% ------ X components ------ %

BOUND_COND = spline(linspace(0,1,nPts),uv(:,1)'); %#ok<NASGU>
%Initialize mesh
[p,e,t] = initmesh('boundaryGeometry');

%Refine and jiggle some before starting adaptive refinement.
for j = 1:min(2,meshQuality)
    [p,e,t] = refinemesh('boundaryGeometry',p,e,t);   
end
p = jigglemesh(p,e,t,'Opt','minimum','Iter',meshQuality*20);

%Do a few rounds of adaptive refinement
[~,px,ex,tx] = adaptmesh('boundaryGeometry','boundaryCondition',...
                         pdePar(1),pdePar(2),pdePar(3),...
                         'Ngen',meshQuality,'Mesh',p,e,t);
%Refine again
for j = 1:min(1,ceil(meshQuality/2))
    [px,ex,tx] = refinemesh('boundaryGeometry',px,ex,tx);   
end
px = jigglemesh(px,ex,tx,'Opt','minimum','Iter',meshQuality*20);

%Get the final solution with a few more rounds of refinement
[uX,px,ex,tx] = adaptmesh('boundaryGeometry','boundaryCondition',...
                           pdePar(1),pdePar(2),pdePar(3),...
                           'Ngen',meshQuality,'Mesh',px,ex,tx);


% ------- Y components --------- %

BOUND_COND = spline(linspace(0,1,nPts),uv(:,2)');

%Do a few rounds of adaptive refinement
[~,py,ey,ty] = adaptmesh('boundaryGeometry','boundaryCondition',...
                        pdePar(1),pdePar(2),pdePar(3),...
                        'Ngen',meshQuality,'Mesh',p,e,t);
%Refine again
for j = 1:min(1,ceil(meshQuality/2))
    [py,ey,ty] = refinemesh('boundaryGeometry',py,ey,ty);   
end
py = jigglemesh(py,ey,ty,'Opt','minimum','Iter',meshQuality*10);

%Get final solution with a few more rounds of adaptive refinement
[uY,py,ey,ty] = adaptmesh('boundaryGeometry','boundaryCondition',...
                          pdePar(1),pdePar(2),pdePar(3),...
                          'Ngen',meshQuality,'Mesh',py,ey,ty);


%% ------- Output ----- %%

if nargout > 2
    %If the solution was requested in mesh form
    varargout{1}=uX;
    varargout{2}=px;
    varargout{3}=ex;
    varargout{4}=tx;
    varargout{5}=uY;
    varargout{6}=py;
    varargout{7}=ey;
    varargout{8}=ty;
else
    %Otherwise, interpolate the mesh solution to gridded data
    varargout{1} = tri2grid(px,tx,uX,1:imgSize(1),1:imgSize(2));
    varargout{2} = tri2grid(py,ty,uY,1:imgSize(1),1:imgSize(2));
end



