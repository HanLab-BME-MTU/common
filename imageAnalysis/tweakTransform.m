function imProd = tweakTransform(baseImage,inImage,dX)

%This function calculates the mean squared error between two images after
%one is transformed by the transformation specifed by the vector dX. The
%vector should be a vectorized version of the transformation matrix as used
%by t imtransform. 

%Assemble the transformation matrix
tMat = zeros(3);

%If limited to homogenous matrices and projective transformations

%Check the number of parameters input and use this to determine the
%transformation type.
if length(dX(:)) == 2
    tMat = eye(3);
    tMat(3,1:2) = dX(:)';
    tType = 'projective';    
elseif length(dX(:)) == 6
    tMat(3,3) = 1;
    tMat(1:6) = dX(:);
    tType = 'projective';
elseif length(dX(:)) == 9;   %If non-homogenous allowed in proj. xform
    tMat(:) = dX(:);
    tType = 'projective';
elseif length(dX(:)) == 20;  %If polynomial transform
    tType = 'polynomial';
    tMat = zeros(10,2);
    tMat(:) = dX(:);
end

%Create the transform
if strcmp(tType,'projective')
    xForm = maketform('projective',tMat);
else
    xForm = cp2tform(ones(10,2),ones(10,2)+rand(10,2)/10,'polynomial',3);
    xForm.tdata = tMat;
end

%Transform the input image with the new transformation
inImage = imtransform(inImage,xForm,'XData',[1 size(baseImage,2)],'YData',[1 size(baseImage,1)],'FillValues',1);

%Calculate the product of the two images. Ignore a border around the edge
%to prevent biasing by the filled-in pixels not present in image
brdr = 20;
imDiff = abs(baseImage(brdr:end-brdr,brdr:end-brdr) - inImage(brdr:end-brdr,brdr:end-brdr));
imProd = sum(sum(imDiff ./ baseImage(brdr:end-brdr,brdr:end-brdr)));
