function fImg=locmax2d(img,mask,keepFlat)
%LOCALMAX searches for local maxima in an image
%
%    SYNOPSIS fImg=(img,mask)
%
%    INPUT    img    image matrix
%             mask   EITHER a vector [m n] that defines the operator window
%                    dimensions (all ones).
%                    OR a structural element, a matrix, that contains only
%                    0/1. Structural elements such as discs can be defined
%                    easily using the built-in matlab function "strel".
%                    The input matrix must have an odd number of columns
%                    and rows. 
%             keepFlat Optional input variable to choose whether to remove
%                      "flat" maxima or to keep them. Default is 0, to
%                      remove them. - KJ
%
%    OUTPUT   fImg   map with all the local maxima (zeros elsewhere);
%                    the non-zero values contain the original value 
%                    of the image at that place
%
% NOTE: convert fImg to uint8 or uint16 for optimal display!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PARAMETER CHECK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
   error('Please define all parameters');
end

% If mask is a vector with two entries, then, these two entries define the
% number of rows and columns of a rectangular mask filled with ones:
if length(mask(:))==2
    % make sure the mask elements are odd numbers (only then, the 
    % local max operator is properly defined)
    indx = find(~mod(mask,2));
    mask(indx) = mask(indx) + 1;
    rows=mask(1);
    cols=mask(2);
    % number of non-zero elements:
    numEl=prod(mask);
    
    % generate the real mask:
    mask=ones(mask);    
else
    % mask is a flat structural element, with an odd number of cols and
    % rows. Note that this excludes the only possible overlap case [1 1] 
    % which would be treated as identical operation above.
    
    % first check that mask has the right size:
    [rows cols]=size(mask);
    if ~mod(rows,2) || ~mod(cols,2)
        % There is no simple way of extending a general mask, thus:
        error('Mask must have an odd number of rows and columns!');
    end
    
    % get matrix dimensions:
    [rows cols]=size(mask);
    
    % Check that the matrix contains only 0/1. This check could be omitted:
    checkMat = (mask==1 | mask==0);
    if sum(checkMat(:))<rows*cols
        % There is no simple way of extending a general mask, thus:
        error('The mask may contain binary values only!');
    end    
    
    % number of non-zero elements:
    numEl=sum(mask(:));
end

if nargin < 3 || isempty(keepFlat)
    keepFlat = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFINITIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply a max filter
fImg = ordfilt2(img,numEl,mask);
if keepFlat == 0 %change made by KJ
    fImg2 = ordfilt2(img,numEl-1,mask);
    fImg(fImg2==fImg)=0;
end

% take only those positions where the max filter and the original image value
% are equal -> this is a local maximum
fImg(fImg ~= img) = 0;

% set image border to zero
auxM = zeros(size(img));
auxM(fix(rows/2)+1:end-fix(rows/2),fix(cols/2)+1:end-fix(cols/2)) = ...
    fImg(fix(rows/2)+1:end-fix(rows/2),fix(cols/2)+1:end-fix(cols/2));
fImg=auxM;
