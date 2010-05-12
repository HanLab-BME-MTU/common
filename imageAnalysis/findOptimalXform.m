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

%Bounds which are common to transformations
xX = [-5 5]; %X translation
yX = [-5 5]; %Y translation
xSc = [.9 1.1]; %Scaling in x direction
ySc = [.9 1.1]; %Scaling in y direction
ySk = [-.1 .1]; %Y "skew"
xSk = [-.1 .1]; %X "skew"
        
switch xType
    
    
    case 'projective'
        %Constraints for projective

        xySc = [.9 1.1]; %Scaling in "z" - x&y simultaneously INVERTED - >1 decreses size!
        
        xL = [-1e-4 1e-4]; %X "lean" - perspective shift + moves right side of image in negative z
        yL = [-1e-4 1e-4]; %Y "lean" - perspective shift + moves bottom side of image in negative z

        Amin = [xSc(1),...
                xSk(1),...
                xX(1),....
                ySk(1),...
                ySc(1),...
                yX(1),...
                xL(1),...
                yL(1),...
                xySc(1)]';

        Amax = [xSc(2),...
                xSk(2),...
                xX(2),....
                ySk(2),...
                ySc(2),...
                yX(2),...
                xL(2),...
                yL(2),...
                xySc(2)]';    
    case 'polynomial'
   
        %Arbitrary small number for higher-order, non-linear terms
        bS = 1e-5;
        
        Amin = [...
            xX(1),yX(1); %X and Y translation            
            xSc(1),ySk(1);%xScaling and y skewing
            xSk(1),ySc(1); %x Skewing and yScaling
            -bS, -bS; %bulging xy and yx
            -bS,  -bS; %barrell x and arc y
            -bS,  -bS; %arc x and barell y
            -5*bS,  -5*bS;    %???? wierd higher order shit...
            -10*bS,  -10*bS;
            -100*bS, -100*bS;
            -1000*bS,-1000*bS;];
        
      Amax = [...
            xX(2),yX(2); %X and Y translation            
            xSc(2),ySk(2);%xScaling and y skewing
            xSk(2),ySc(2); %x Skewing and yScaling
            bS, bS; %bulging xy and yx
            bS,  bS; %barrell x and arc y
            bS,  bS; %arc x and barell y
            5*bS,  5*bS;    %???? wierd higher-order shit...
            10*bS,  10*bS;
            100*bS,  100*bS;
            1000*bS, 1000*bS;];
             
        
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
    
    baseImage(baseImage == 0) = min(baseImage(baseImage>0));
    inputImage(inputImage == 0) = min(inputImage(inputImage>0));
    
end



%Create objective function for minimization
objFun = @(x)tweakTransform(baseImage,inputImage,x);

%minOptions = optimset('MaxFunEvals',5e3);
tic;
%Minimize the objective function to find optimal transform
disp('Please wait, calculating transform...');

minOpts = optimset('TolFun',1e-8);


[x,fval,exflag,output] = fmincon(objFun,iGuess(:),[],[],[],[],Amin,Amax,[],minOpts);
%x = fminsearch(objFun,iGuess(:));
%x = fminbnd(objFun,Amin,Amax);

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
    
    
    fsFigure(.75);
    subplot(1,2,1)
    image(cat(3,mat2gray(baseImage),mat2gray(inputImage),zeros(size(inputImage))))
    title('Original Overlay'),axis image
    newImage = imtransform(inputImage,optimalXform,'XData',[1 size(baseImage,2)],'YData',[1 size(baseImage,1)],'FillValues',1);
    subplot(1,2,2)
    image(cat(3,mat2gray(baseImage),mat2gray(newImage),zeros(size(inputImage))))
    title('Aligned Overlay'),axis image
    
end

function imErr = tweakTransform(baseImage,inImage,dX)

%This function calculates the mean squared error between two images after
%one is transformed by the transformation specifed by the vector dX. The
%vector should be a vectorized version of the transformation matrix as used
%by imtransform.m 

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
if strcmp(tType,'projective') %THIS IS WASTEFUL! set tForm as global variable
    xForm = maketform('projective',tMat);
else
    xForm = cp2tform(ones(10,2),ones(10,2)+rand(10,2)/10,'polynomial',3); 
    xForm.tdata = tMat;
end

%Transform the input image with the new transformation
inImage = imtransform(inImage,xForm,'XData',[1 size(baseImage,2)],'YData',[1 size(baseImage,1)],'FillValues',NaN);

%Calculate the error between the two images. Ignore NaNs.
%imErr = sqrt(nanmean(((baseImage(:) - inImage(:)) .^2 ) ./ baseImage(:) ));
imErr = sqrt(nanmean(((baseImage(:) - inImage(:)) .^2 )));
