function dataM = imDataMap(imgDim,x,y,data,varargin)
%dataM: Create an image speed map (pixel wise) from sampled speed data.
%
% SYNOPSIS : 
%    dataM = imSpeedMap(imgDim,x,y,data)
%    dataM = imSpeedMap(imgDim,x,y,data,bnd)
%    dataM = imSpeedMap(imgDim,x,y,data,bnd,method)
% INPUT :
%    imgDim : The dimension of the image in the form of [m,n] where m and n is 
%             the vertical and horizontal dimension respectively.
%    x, y   : The x and y coordinates of sampled points.
%    data   : The speed at sampled points.
%
%    Optional inputs:
%    bnd    : An Nx2 matrix that specifies a polygon bondary for the sampling
%             points. 'bnd(:,1)' are the x-coordinates and 'bnd(:,2)' are the
%             y-coordinates. We only calculate the speed by interpolation for 
%             pixels inside this polygon. The default is to calculate for all 
%             pixels.
%    method : String that specifies the interpolation method:
%       'linear'    - Triangle-based linear interpolation (default).
%       'cubic'     - Triangle-based cubic interpolation.
%       'nearest'   - Nearest neighbor interpolation.
%       'v4'        - MATLAB 4 griddata method.
%       The default is 'nearest'.
%
% OUTPUT :
%    dataM : A pixel wise map of the data of dimension 'imgDim'.

%The defaults.
bnd    = [];
method = 'nearest';

if nargin > 6
   error('Too many input arguments.');
elseif nargin > 4
   bnd = varargin{1};
elseif nargin > 5
   method = varargin{2};
end

%Create grids for pixels.
[gridX,gridY] = meshgrid(1:imgDim(2),1:imgDim(1));
dataM = griddata(x,y,data,gridX,gridY,method);

%Set speed at pixels out of 'bnd' to be NaN.
if ~isempty(bnd)
   [in on] = inpolygon(gridX(:),gridY(:),bnd(:,1),bnd(:,2));
   outI = find(in==0 | on==1);
   dataM(outI) = NaN;
end

