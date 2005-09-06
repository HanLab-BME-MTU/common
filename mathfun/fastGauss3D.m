function out=fastGauss3D(img,sigma,fSze,correctBorder,filterMask);
% fastGauss3D	apply a 3 dimensional gauss filter
%
%    SYNOPSIS out=Gauss3D(img,sigma,fSze);
%
%    INPUT:   img      3-dimensional data
%             sigma    of gauss filter [sigmaX sigmaY sigmaZ]
%             fSze     size of the gauss mask [sizeX sizeY sizeZ]
%                          (odd size required for symmetric mask!)
%             correctBorder (optional) if 1, border effects of the filtering are
%                              lessened
%             filterMask    (optional) supply your own mask instead of a
%                              Gauss. In this case, sigma and fSze don't
%                              have to be supplied
%
%    OUTPUT:  out      filtered data

% c: 13/03/01 dT

if nargin > 4 && ~isempty(filterMask)
    % Use user-supplied filter
    gauss = filterMask;
    sigma = [];
    fSze = size(filterMask);
    filterMask = []; % reclaim memory
else
    % make a GaussFilter
    gauss=GaussMask3D(sigma,fSze);
end
    

convnOpt = 'same';

if nargin > 3 & ~isempty(correctBorder)
    if correctBorder == 1
        %correct for border effects
        try
            img = addBorder(img,fSze);
            convnOpt = 'valid';
        catch
            err = lasterr;
            if findstr(err,'image too small')
                disp(['Warning: ',err,' Filtering is likely to produce useless results!'])
                convnOpt = 'same'; %this will work, but is not likely to succeed
            else
                rethrow(lasterr);
            end
        end
    end
end

% Convolve matrices
out=convn(img,gauss,convnOpt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newImg = addBorder(img,fSze);
%adds a border of halFsze around the img, and fill it with pixels whose
%value is computed as follows:
%for all 3 dimensions i (=6 sides) of the array, the median of the first or last
%halFsze(i)+1 pixels in this dimension, respectively, is computed. The
%edges and corners are filled with means of the adjacent sides or edges,
%respectively
%c: 030514 jonas

%init
halFsze = zeros(1,3);
halFsze(1:length(fSze)) = floor(fSze/2); %half filter size
oldSize = size(img);
if length(oldSize) < 3 
    error('image too small for 3D filtering: no 3rd dimension!')
elseif any(oldSize < halFsze+1)
    error('image too small for 3D filtering: filter is larger than image!')
end

newImg = zeros(oldSize+2*halFsze);
newImg(halFsze(1)+1:end-halFsze(1),halFsze(2)+1:end-halFsze(2),halFsze(3)+1:end-halFsze(3)) = img;

%alternatively, this could be made with feval and for-loop, but it would be a
%chore, too

%2 sides first direction
inSide11 = img(1:halFsze(1)+1,:,:);
side11val = median(inSide11(:));
newImg(1:halFsze(1),halFsze(2)+1:end-halFsze(2),halFsze(3)+1:end-halFsze(3)) = side11val;

inSide12 = img(end-halFsze(1):end,:,:);
side12val = median(inSide12(:));
newImg(end-halFsze(1)+1:end,halFsze(2)+1:end-halFsze(2),halFsze(3)+1:end-halFsze(3)) = side12val;

%2 sides second direction
inSide21 = img(:,1:halFsze(2)+1,:);
side21val = median(inSide21(:));
newImg(halFsze(1)+1:end-halFsze(1),1:halFsze(2),halFsze(3)+1:end-halFsze(3)) = side21val;

inSide22 = img(:,end-halFsze(2):end,:);
side22val = median(inSide22(:));
newImg(halFsze(1)+1:end-halFsze(1),end-halFsze(2)+1:end,halFsze(3)+1:end-halFsze(3)) = side22val;

%2 sides third direction
inSide31 = img(:,:,1:halFsze(3)+1);
side31val = median(inSide31(:));
newImg(halFsze(1)+1:end-halFsze(1),halFsze(2)+1:end-halFsze(2),1:halFsze(3)) = side31val;

inSide32 = img(:,:,end-halFsze(3):end);
side32val = median(inSide32(:));
newImg(halFsze(1)+1:end-halFsze(1),halFsze(2)+1:end-halFsze(2),end-halFsze(3)+1:end) = side32val;

%4 edges parallel to first direction
edge2131val = mean([side21val,side31val]);
newImg(:,1:halFsze(2),1:halFsze(3)) = edge2131val;

edge2231val = mean([side22val,side31val]);
newImg(:,end-halFsze(2)+1:end,1:halFsze(3)) = edge2231val;

edge2132val = mean([side21val,side32val]);
newImg(:,1:halFsze(2),end-halFsze(3)+1:end) = edge2132val;

edge2232val = mean([side22val,side32val]);
newImg(:,end-halFsze(2)+1:end,end-halFsze(3)+1:end) = edge2232val;

% 4 edges parallel to second direction
edge1131val = mean([side11val,side31val]);
newImg(1:halFsze(1),:,1:halFsze(3)) = edge1131val;

edge1231val = mean([side12val,side31val]);
newImg(end-halFsze(1)+1:end,:,1:halFsze(3)) = edge1231val;

edge1132val = mean([side11val,side32val]);
newImg(1:halFsze(1),:,end-halFsze(3)+1:end) = edge1132val;

edge1232val = mean([side12val,side32val]);
newImg(end-halFsze(1)+1:end,:,end-halFsze(3)+1:end) = edge1232val;

% 4 edges parallel to third direction
edge1121val = mean([side11val,side21val]);
newImg(1:halFsze(1),1:halFsze(2),:) = edge1121val;

edge1221val = mean([side12val,side21val]);
newImg(end-halFsze(1)+1:end,1:halFsze(2),:) = edge1221val;

edge1122val = mean([side11val,side22val]);
newImg(1:halFsze(1),end-halFsze(2)+1:end,:) = edge1122val;

edge1222val = mean([side12val,side22val]);
newImg(end-halFsze(1)+1:end,end-halFsze(2)+1:end,:) = edge1222val;

%corner 000
newImg(1:halFsze(1),1:halFsze(2),1:halFsze(3)) = mean([edge2131val,edge1131val,edge1121val]);

%corner 100
newImg(end-halFsze(1)+1:end,1:halFsze(2),1:halFsze(3)) = mean([edge2131val,edge1231val,edge1221val]);

%corner 010
newImg(1:halFsze(1),end-halFsze(2)+1:end,1:halFsze(3)) = mean([edge2231val,edge1131val,edge1122val]);

%corner 110
newImg(end-halFsze(1)+1:end,end-halFsze(2)+1:end,1:halFsze(3)) = mean([edge2231val,edge1231val,edge1222val]);

%corner 001
newImg(1:halFsze(1),1:halFsze(2),end-halFsze(3)+1:end) = mean([edge2132val,edge1132val,edge1121val]);

%corner 101
newImg(end-halFsze(1)+1:end,1:halFsze(2),end-halFsze(3)+1:end) = mean([edge2132val,edge1232val,edge1221val]);

%corner 011
newImg(1:halFsze(1),end-halFsze(2)+1:end,end-halFsze(3)+1:end) = mean([edge2232val,edge1132val,edge1122val]);

%corner 111
newImg(end-halFsze(1)+1:end,end-halFsze(2)+1:end,end-halFsze(3)+1:end) = mean([edge2232val,edge1232val,edge1222val]);