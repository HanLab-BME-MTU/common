function skel = skeleton(bw)
%SKELETON Morphological binary skeletonization in 2D or 3D
% 
% skel = skeleton(BW)
% 
% This function performs skeletonization on the input (2D or 3D) binary
% matrix using erosion and opening [1] p 387. PLEASE NOTE: Due to the
% inherent discretization, this method of skeletonization DOES NOT
% GUARANTEE connectedness of the resulting skeleton [1] p. 386 Fig. 9.33,
% even if the original object is connected.
% 
%
% Reference: 
% [1] A.K. Jain - Fundamentals of Digital Image Processing, Prentice-Hall,
% 1989
%
% Input:
% 
%   bw - The 2D or 3D binary matrix to skeletonize.
% 
% 
% Output:
% 
%   skel - The 2D or 3D binary matrix with the skeleton points.
% 
% Hunter Elliott 
% 2/2010
%

if ndims(bw) == 3
    nH_1 = binarySphere(1);
elseif ndims(bw) == 2
    nH_1 = strel('disk',1,0).getnhood;
else
    error('Input matrix must be 2D or 3D!!!')
end

%Make sure the matrix is logical
bw = logical(bw);

%Initialize the skeleton matrix
skel = false(size(bw));

%Initialze the eroded matrix
imE = true(size(bw));

j = 1;
while any(imE(:))
    
    if ndims(bw) == 3
        nH = binarySphere(j);
    else
        nH = strel('disk',double(j),0).getnhood;
    end    
    %perform the erosion at the current radius.
    imE = imerode(bw,nH);
    %Subtract the opening with unit radius.
    skel = skel | imE - imopen(imE,nH_1);
    
    j = j + 1;
    
end

