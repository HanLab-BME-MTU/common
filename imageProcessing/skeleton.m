function skel = skeleton(bw,method,divThreshold)
%SKELETON Morphological binary skeletonization in 2D or 3D
% 
% skel = skeleton(BW)
% 
% skel = skeleton(BW,'e')
%
% skel = skeleton(BW,'d',divThreshold);
%
% This function performs skeletonization on the input (2D or 3D) binary
% matrix, and supports a variety of methods for performing this operation.
%
%
% Input:
% 
%   bw - The 2D or 3D binary matrix to skeletonize.
% 
%   method - Character array with the name of the method to use for
%   skeletonization one of these possibilities:
% 
%       'e' - Erosion & opening. As described in [1] p387. PLEASE NOTE: Due to
%       the inherent discretization, this method of skeletonization DOES
%       NOT GUARANTEE connectedness of the resulting skeleton [1] p. 386
%       Fig. 9.33, even if the original object is connected. However, it IS
%       guaranteed that the returned points will be on the medial axis /
%       medial surface of the input mask.
% 
%       't' - Thinning. The skeleton is produced by repeated thinning of
%       the input binary image until no more points are removed. For 2D
%       objects this produces results very similar to bwmorph('skel'),
%       while for 3D objects it produces a reasonable approximation of the
%       medial surface. Connectivity is preserved. This approach is very
%       slow!
% 
%       'd' - Divergence. This method uses the divergence of the gradient
%       of the distance transform of the input binary image to find
%       skeleton points. Connectivity is NOT guaranteed to be preserved,
%       and the resulting points are not guaranteed to be on the true
%       medial axis. HOWEVER, this method is fast and is much less
%       sensitive to minor variations in the shape of the input binary
%       object than the methods above. The skeleton lies on ridge-lines of
%       the distance transform. These ridge lines are areas where the
%       divergence of the gradient of the distance transform is negative.
%       The skeleton is returned by thresholding this divergence at the
%       value divThreshold:
%
%           divThreshold - The value to threshold the divergence at to
%           produce the skeleton (should be negative!). Lower (more
%           negative) values will give fewer points which are more likely
%           to lie on the true keleton, but less likely to preserve
%           connectivity. Optional. Default is -1;
%
%       Optional. Default method is 't'
%
%
% Output:
% 
%   skel - The 2D or 3D binary matrix with the skeleton points.
% 
%
% References: 
% [1] A.K. Jain - "Fundamentals of Digital Image Processing", Prentice-Hall,
% 1989
%
%
% Hunter Elliott 
% 2/2010
%

%% ------ Input ----- %%
if nargin < 1 || isempty(bw) || ndims(bw) < 2 || ndims(bw) > 3
    error('Must input a 2D or 3D matrix for skeletonization!')
end

if nargin < 2 || isempty(method)
    method = 't';
end

if nargin < 3 || isempty(divThreshold)
    divThreshold = -1;
elseif divThreshold >= 0
    warning('Non-negative threshold values are not recommended! What you''re getting is NOT a skeleton!') %#ok<WNTAG>
end
    

%Make sure the matrix is logical
bw = logical(bw);

switch method
    
    
    case 'e'
    
        % ---- Skeletonization by Erosion & Opening ---- %
        
        if ndims(bw) == 3
            nH_1 = binarySphere(1);
        else
            nH_1 = strel('disk',1,0).getnhood;
        end
        
        %Initialize the skeleton matrix
        skel = false(size(bw));

        j = 1;
        while any(bw(:))

            if ndims(bw) == 3
                nH = binarySphere(j);
            else
                nH = strel('disk',double(j),0).getnhood;
            end    
            %perform the erosion at the current radius.
            bw = imerode(bw,nH);
            %Subtract the opening with unit radius.
            skel = skel | bw ~= imopen(bw,nH_1);

            j = j + 1;
            
        end
        
        
    case 't'
        
        % ----- Skeletonization by Thinning ----- %
        
        %Determine number of non-zero pixels/voxels
        nnzNew = nnz(bw);
        nnzPrev = Inf;
        skel = bw;
        
        %Thin until no more pixels are removed
        while nnzNew < nnzPrev
            
            nnzPrev = nnzNew;
            
            skel = bwthin(skel);
        
            nnzNew = nnz(skel);
            
        end
        
        
    case 'd'
        
        % ----- Skeletonization by Divergence ----- %
        
        %Calculate the gradient of the distance transform of the inverse of
        %the input object
        [gX,gY,gZ] = gradient(bwdist(~bw));
        
        %Calculate the divergence of this gradient vector field and
        %threshold it.
        skel = divergence(gX,gY,gZ) < divThreshold;                
        
        
    otherwise
        
        error(['The input "' method '" is not a recognized method!'])
        
end


