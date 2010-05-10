function optimalXform = findOptimalXform(baseImage,inputImage,showFigs,xType)
%FINDOPTIMALXFORM finds a transform which maximizes the agreement between the input images
%
% optimalXform = findOptimalXform(baseImage,inputImage,showFigs,xType)
%
% This function finds a transform which maximizes the agreement between the
% input and base images by transforming the input image. This transform can
% then be used to register two channels in a movie.
% 
% Input:
% 
%  baseImage - The image the input image will be compared (registered) to.
% 
%  inputImage - The image to transform to match the base image.
%
%  showFigs - True or false. If true, figures will be show.
%
%  xType - A characted string describing the transformation type to use to
%  align the two images. Default is 'projective'
% 
%       Available Transformation Types:
%
%           'projective' - Projective transformation.
%
%           'projectiveH' - Homogeneous projective transformation
%
%           'polynomial' - Non-linear, polynomial transform
%
%           'translation' - X-Y translation only.
%
%Hunter Elliott, 2008
%

if nargin < 3 || isempty(showFigs)
    showFigs = 0;
end

if nargin < 4 || isempty(xType)
    xType = 'projective';
end



maxDisp = 20; %constraint for displacement - maximum number of pixels to displace

switch xType

    case 'projective'    
        iGuess = eye(3);        
    case 'projectiveH' %homogeneous projective transform
        iGuess = [1 0 0 0 1 0];                
    case 'polynomial'
        iGuess = zeros(10,2);
        iGuess(2) = 1;
        iGuess(13) = 1;
    case 'translation'
        iGuess = [0 0];
end

%If the images are actually masks, scale them differently to avoid rounding
%effects
if islogical(baseImage) && islogical(inputImage)
    baseImage = cast(baseImage,'double');
    inputImage = cast(inputImage,'double');
    baseImage = baseImage ./ 2 + .1;
    inputImage = inputImage ./ 2 + .1;    
else
    %Normalize the images
    baseImage = cast(baseImage,'double');
    inputImage = cast(inputImage,'double');
    baseImage = baseImage - min(baseImage(:));
    baseImage = baseImage ./ max(baseImage(:));
    inputImage = inputImage - min(inputImage(:));
    inputImage = inputImage ./ max(inputImage(:));
end



%Create objective function for minimization
objFun = @(x)tweakTransform(baseImage,inputImage,x);

minOptions = optimset('MaxFunEvals',5e3);
tic;
%Minimize the objective function to find optimal transform
disp('Please wait, calculating transform...');
[x,fval,exflag,output] = fmincon(objFun,iGuess(:),ones(1,length(iGuess(:))),maxDisp,[],[],[],[],[],minOptions);
%x = fminsearch(objFun,iGuess(:));
telaps = toc;
if exflag > 0
    disp(['Finished. Took ' num2str(telaps/60) ' minutes.']);
else
    disp('Optimization failed!')
    optimalXform = [];
    return
end


if strcmp(xType,'translation')
    tMat = eye(3);
    tMat(3,1:2) = x(:)';    
else
    tMat = zeros(size(iGuess));
    tMat(:) = x(:);        
end
if strcmp(xType,'projective') || strcmp(xType,'translation')
    optimalXform = maketform('projective',tMat);
else
    optimalXform = cp2tform(ones(10,2),ones(10,2)+rand(10,2)/10,'polynomial',3);
    optimalXform.tdata = tMat;
end

if showFigs
    
    
    scrSize = get(0,'ScreenSize');
    resFig = figure('position',[75 75 round(scrSize(3)-200) scrSize(4)-200]);
    
    subplot(1,2,1)
    imagesc(baseImage ./ inputImage)
    title('Original Ratio'),axis image,caxis([.5 2]);
    newImage = imtransform(inputImage,optimalXform,'XData',[1 size(baseImage,2)],'YData',[1 size(baseImage,1)],'FillValues',1);
    subplot(1,2,2)
    imagesc(baseImage ./ newImage)
    title('Aligned Ratio'),axis image,caxis([.5 2]);
    
end

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
