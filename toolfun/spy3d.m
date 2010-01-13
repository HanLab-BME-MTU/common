function spy3d(matIn)

%
% spy3d(matIn)
%
% The spy function, but in 3d!!
% (This function visualizes the sparseness of an input 3d matrix.)
%
% Input:
%
% matIn - The matrix to visualize sparseness of. Must be 3D.
%
%
% Output:
%
%   No output, just a 3D figure with dots in each voxel which has a
%   non-zero value. The color of these dots transitions from red to blue
%   in the z-direction.
%
%
% Hunter Elliott
% Sometime in 2008?
%


if ndims(matIn) ~= 3
    error('Input matrix must be 3d!!')
end
    

P = size(matIn,3);

hold on

for p = 1:P
    
    %Find the non-zero values
    [row,col] = find(matIn(:,:,p));
    
    %Plot them at the appropriate z-location
    plot3(col,row,p * ones(1,length(row)),'.','color',[ (P-p)/P 0 p/P ]  );
    
    
end

%Set the viewpoint to 3d
view(3)
%Set axes scaling
axis image