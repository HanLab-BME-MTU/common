function dataM = imDataMap(imgDim,YX,data,varargin)
%imDataMap: Create a pixelwise data map from sampled data on either grid or 
%           irregular data points.
%
% One problem with interpolation of sampled data values to pixels of an
% underlying image is that the data set does not cover the whole image.
% In this case, one wants to distinguish between the influence region of
% the sampled data and the region without data. Here, we use 'bwdist' to
% first calculate the nearest distance matrix of each pixel to the sampled
% data points. Then, we set pixels whose nearest distances are greater than
% an influence length to be no data pixels.
% 
% SYNOPSIS : 
%    dataM = imDataMap(imgDim,YX,data)
%    dataM = imDataMap(imgDim,YX,data,'par1',value1,...)
%
% INPUT :
%    imgDim : The dimension of the image in the form of [m,n] where m and n is 
%             the vertical and horizontal dimension respectively.
%    YX     : The y and x coordinates of sampled points. If it is a cell of
%             two vectors in the form {y,x}, a mesh grid will be generated.
%             Otherwise, it is expected to be a two-column matrix in the
%             form [y,x] that gives the coordinate of each data point.
%    data   : The sampled data. When YX = {y,x} is a cell, it is a 2D matrix of
%             size ny-by-nx where ny is the length of y and nx is the 
%             length of x. When YX is a two-column matrix, it is a 1D vector
%             of data values.
%
%    Optional parameter pairs:
%    infLen : The influence length that specifies the radius of a disk
%             around each sampled data point. The union of the disks defines 
%             the influence region of the data points. Pixels outside the 
%             influence region are considered no data points which are
%             represented by 'NaN' in the output 'dataM'. By default, it is
%             the average of the grid intervals in both x andy direction if
%             the data points are grids, or it is the average of the
%             nearest distances of each data point if the data points are
%             irregular.
%    bnd    : An Nx2 matrix that specifies a polygon bondary for the sampling
%             points. 'bnd(:,1)' are the x-coordinates and 'bnd(:,2)' are the
%             y-coordinates. We only interpolate for pixels inside this polygon. 
%             The default is to calculate for all pixels.
%    mask   : A black-white mask can also be given. Then, only pixels with
%             values of 1 in the mask are interpolated. It has to be of
%             size 'imgDim'.
%    method : String that specifies the interpolation method for griddata:
%       'linear'    - Triangle-based linear interpolation (default).
%       'cubic'     - Triangle-based cubic interpolation.
%       'nearest'   - Nearest neighbor interpolation.
%       'v4'        - MATLAB 4 griddata method.
%       The default is 'nearest'.
%    order  : The order of the spline when the data are on grids.
%
% OUTPUT :
%    dataM : A pixel wise map of the data of dimension 'imgDim'.
%
% By Lin Ji.

%Get the dimension of the image.
m = imgDim(1);
n = imgDim(2);

%The defaults.
infLen = [];
bnd    = [];
mask   = [];
method = 'nearest';
order  = 4;

if iscell(YX)
    if length(YX) > 2 | ~isnumeric(YX{1}) | ~isnumeric(YX{2})
        error(['The coordinates of data sampled on grids ' ...
            'should be entered in the form {y,x}.']);
    end
    y = YX{1};
    x = YX{2};
    
    if ~isempty(find(diff(y)<=0)) | ~isempty(find(diff(x)<=0))
        error('The grids in each dimension have to be increasing sequence.');
    end
    
    if size(data,1) ~= length(y) | size(data,2) ~= length(x)
        error('The size of the data does not match the size of the grids.');
    end
    
    %Exclude grids out of image area.
    y(find(y<1 | y>m)) = [];
    x(find(x<1 | x>n)) = [];
    
    if isempty(y) | isempty(x)
        error('No data points are inside the image area.');
    end
    
    data = data(find(y>=1 & y<=m),find(x>=1 & x<= n));
    
    %Data is on grids.
    dataLoc = 'grid';
    
    
    %Create grid points.
    [X,Y] = meshgrid(x,y);
elseif isnumeric(YX)
    if ndims(YX) > 2 | size(YX,2) ~= 2
        error(['The coordinates of randomly sampled data sets ' ...
            'should be entered as a two-column matrix.']);
    end
    
    %Exclude data points are all outside the image area
    outInd = find(YX(:,1)<1 | YX(:,1)>m | YX(:,2)<1 | YX(:,2)>n);
    YX(outInd,:) = [];
    data(outInd) = [];
    
    if isempty(YX)
        error('No data points are inside the image area.');
    end
    
    Y = YX(:,1);
    X = YX(:,2);
    
    %Data locations are random (irregular).
    dataLoc = 'random';
else
    error('The coordinates of the sampled data are not recogonized.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parsing optional parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 3
   if mod(nargin-3,2) ~= 0
       error('Wrong number of arguments. Parameter/value has to be in pair.');
   end
   for k = 1:2:nargin-3
       switch varargin{k}
           case 'infLen'
               infLen = varargin{k+1};
           case 'bnd'
               bnd = varargin{k+1};
           case 'mask'
               mask = varargin{k+1};
           case 'method'
               method = varargin{k+1};
           case 'order'
               order = varargin{k+1};
       end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the infulence region of the data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numInd = find(~isnan(data));
if isempty(numInd)
    dataM = NaN*ones(m,n);
    return;
end

if isempty(infLen)
    %The user does not specify 'infLen'. We set the default length to be
    % the mean distance between sampled data points.
    if strcmp(dataLoc,'grid')
        infLen = (mean(diff(y))+mean(diff(x)))/2;
    else
        %When the data points are randomly distributed, we use the average
        % of the nearest distance for each point as the default 'infLen'.
        % To do this, we need to calculate the distance matrix between
        % points.
        eps = 1e-10; %See help creatSparseDistanceMatrix.
        D = createSparseDistanceMatrix(YX(numInd,:),YX(numInd,:),m+n,eps);
        infLen = 0;
        for k = 1:size(D,1)
            infLen = infLen+min(D(k,find(D(k,:)>eps)),[],2);
        end
        infLen = infLen/size(D,1);
    end
end

%In order to use 'bwdist', we first create a mask where only sampled data
% pixels are nonzero.
bwData = zeros(m,n);
bwData((round(X(numInd))-1)*m+round(Y(numInd))) = 1;

%Calculate the 2D distance transform.
DT = bwdist(bwData,'chessboard');

%Index of pixels that are in the influnce region of the data set.
interpPixelsInd = find(DT<=infLen);

%Create a mask of the influnce region where data are interpolated. We keep
% the variable name 'bwData'.
bwData(:) = 0;
bwData(interpPixelsInd) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(dataLoc,'grid')
    %In this case, we use splines (spapi). Before we use splines, we need
    % first to assign some data value to grids that do not have sampled
    % data. For this we use 'nearest' neighbor interpolation. This is just
    % to guanrantee a smooth extropolation at the boundary between the data
    % region and no-data region.
    nanInd = find(isnan(data));
    numInd = find(~isnan(data));
    
    if ~isempty(nanInd)
        data(nanInd) = griddata(Y(numInd),X(numInd),data(numInd), ...
            Y(nanInd),X(nanInd),'nearest');
    end
    
    %Spline interpolation.
    sp = spapi({order,order},{y,x},data);
    
    %Calculate the data map.
    dataM = fnval(sp,{1:m,1:n});
else
    %In this case we use 'griddata' with the specified 'method'.
    [pixelX,pixelY] = meshgrid(1:n,1:m);
    dataM = griddata(Y,X,data,pixelY,pixelX,method);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use mask 'bwData' to cut off no-data region.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataM(find(bwData==0)) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use user entered mask 'mask' to cut off region of no interest.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(mask)
    dataM(find(mask==0)) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set dataat pixels out of user defined 'bnd' to be NaN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(bnd)
   [in on] = inpolygon(gridX(:),gridY(:),bnd(:,1),bnd(:,2));
   outI = find(in==0 | on==1);
   dataM(outI) = NaN;
end

