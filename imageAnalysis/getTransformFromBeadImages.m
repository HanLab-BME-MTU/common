function xForm = transformFromBeadImages(baseImage,inputImage,xFormType,beadRad,showFigs)
%TRANSFORMFROMBEADIMAGES calculates an alignment transform based on two bead images
% 
% xForm = transformFromBeadImages(baseImage,inputImage)
% xForm = transformFromBeadImages(baseImage,inputImage,xFormType)
% xForm = transformFromBeadImages(baseImage,inputImage,xFormType,beadRad)
% 
% This function is used to calculate a transform which can be used to align
% images from two different cameras/channels on the same microscope. It
% takes as input images of fluorescent beads taken from the two different
% channels and returns a transform which can be used with imtransform.m to
% align the two channels. This is accomplished by detecting bead locations
% in each image and then using cp2tform.m
% 
% ***NOTE*** This function expects that the beads will be sparse within the
% image. That is, the vast majority of the image will be background, and
% the maximum misalignment between the two images is much less than the
% average spacing between beads (at least a factor of ten)
%
% Input:
% 
%   baseImage - The 2D image which is used as a reference for transformation.
%   That is, the transformation will be used to align images to this image.
% 
%   inputImage - The 2D image which will be aligned to the baseImage by the
%   transformation. When the resulting transform is applied to this image,
%   it should align with the baseImage.
% 
%   xFormType - Character array. Optional. Specifies the type of transform
%   to use to align the two images. Default is 'projective', but any type
%   supported by imtransform.m can be used.
%
%   beadRad - Scalar. The approximate radius of the beads in each image,
%   in pixels. Optional, default is 3 pixels.
%
%   showFigs - True/False. If true, figures showing the bead detection and
%   alignment willb e displayed. Optional. Default is false.
%
% Output:
%
%   xForm - The structure describing the transform, as used by
%   imtransform.m
%
% Hunter Elliott 
% 10/2010
%
%% -------- Parameters ------------ %%

spacingFactor = 10; %The spacing between beads in a single image must be at least this factor larger than the maximum bead mis-alignment between images.


%% ------------ Input ----------- %%

if nargin < 2 || isempty(baseImage) || isempty(inputImage)
    error('You must input a base and input image!')
end

[M,N] = size(baseImage);
if size(inputImage,1) ~= M || size(inputImage,2) ~= N
    error('The base and input images must be the same size!');
end

if nargin <3 || isempty(xFormType)
    xFormType = 'projective';
end

if nargin < 4 || isempty(beadRad)
    beadRad = 3;
end

if nargin < 4 || isempty(showFigs)
    showFigs = false;
end

%% ---------- Init -------------- %%

%Determine filter size for local maxima detection
fSize = round(beadRad*2+1);


%% ---------- Detection --------- %%

disp('Detecting beads in both images...')

%Detect local maxima in both images
bMax = locmax2d(baseImage,[fSize fSize]);
iMax = locmax2d(inputImage,[fSize fSize]);

%Keep only the "bright" local maxima. The image should be mostly
%background, so we just use the average of the whole image to get
%approximate background intensity and standard deviation.
bMax = bMax > (mean(baseImage(:))+2*std(double(baseImage(:))));
iMax = iMax > (mean(inputImage(:))+2*std(double(inputImage(:))));





%% ------------ Assignment --------- %%
%Finds correspondence between detcted beads in each image

disp('Determining bead correspondence between images...')

[xB,yB] = find(bMax);
[xI,yI] = find(iMax);

%Number of spots detected in base image
nDetBase = numel(xB);
%Anonymous function for distance calcs
dFun = @(x)(sqrt((x(1)-xI).^2 + (x(2)-yI).^2));

minD = zeros(nDetBase,1);
avgD = zeros(nDetBase,1);
iClosest = zeros(nDetBase,1);

for j = 1:nDetBase

    %Calculate distance from this point to all those in input image
    currDists = dFun([xB(j) yB(j)]);
    
    %find the closest point in the input image and it's distance
    [minD(j),iClosest(j)] = min(currDists);
    
    %Find Average distance
    avgD(j) = mean(currDists);
            
        
end

%TEMP - compare average minimum distance and average spacing to determine quality of detection/image?- HLE

%Get average spacing of all points
avgD = mean(avgD);

%Remove points which do not have a close detection in the input image
xB(minD>(avgD/spacingFactor)) = [];
iClosest(minD>(avgD/spacingFactor)) = [];
yB(minD>(avgD/spacingFactor)) = [];

%Make sure two beads in the second image werent assigned to the same bead
%in the base image
if numel(unique(iClosest)) ~= numel(iClosest)
    error('Problem assigning beed correspondence! Beads may be too closely spaced relative to image misalignment!')
end

%Get ordered list of corresponding points for remaining detections
xI = xI(iClosest);
yI = yI(iClosest);

%Now do sub-resulution detection for the selected points in each image
disp('Performing sub-resolution position refinement...');

%Get window size for gaussian fitting.
winSize = ceil(beadRad)+2;
%Make sure none of the points are too close to the image border
isOK = (xI - winSize) > 0 & (yI - winSize) > 0 & ...
       (xI + winSize) <= M & (yI + winSize) <= M & ...
       (xB - winSize) > 0 & (yB - winSize) > 0 & ...
       (xB + winSize) <= M & (yB + winSize) <= M;
xI = xI(isOK);
yI = yI(isOK);
xB = xB(isOK);
yB = yB(isOK);

nBeads = numel(xB);
       
for j = 1:nBeads
        
    %Get the sub-region of the image for readability/debugging    
    imROI = double(inputImage(xI(j)-winSize:xI(j)+winSize,...
                              yI(j)-winSize:yI(j)+winSize));
    
    %Fit a gaussian to get the sub-pixel location of each bead.
    pVec = fitGaussian2D(imROI,...
                    [0 0 double(inputImage(xI(j),yI(j))) ...
                    beadRad/2 mean(double(inputImage(:)))],'xyAsc');
        
    %Make sure the fit converged, and then add it to the integer position
    if all(abs(pVec(1:2))) < winSize;                  
        xI(j) = xI(j) + pVec(2);
        yI(j) = yI(j) + pVec(1);
    end
    
    %Get the sub-region of the image for readability/debugging  
    imROI = double(baseImage(xB(j)-winSize:xB(j)+winSize,...
                              yB(j)-winSize:yB(j)+winSize));

    %Fit a gaussian to get the sub-pixel location of each bead.                      
    pVec = fitGaussian2D(imROI,...
                [0 0 double(baseImage(xB(j),yB(j))) ...
                 beadRad/2 mean(double(baseImage(:)))],'xyAsc');
        
    if all(abs(pVec(1:2))) < winSize;                             
        xB(j) = xB(j) + pVec(2);
        yB(j) = yB(j) + pVec(1);
    end
     
end



%% ----------- Get Transform ---------- %%

disp('Determining transformation...')

%Use matlab built in function to determine transform
xForm = cp2tform([yI xI],[yB xB],xFormType);


if showFigs
       
    figure
    imshow(cat(3,mat2gray(baseImage),mat2gray(inputImage),zeros(size(baseImage))),[]);
    hold on
    spy(bMax,'ro',beadRad*2+2)
    spy(iMax,'go',beadRad*2+2)       
    
    for j = 1:nBeads
       
        plot([yI(j) yB(j)],[xI(j) xB(j)],'--b')                
        
    end
    
    legend('Base Image','Input Image','Correspondence');    
    
    figure
    xIn = imtransform(inputImage,xForm,'XData',[1 size(baseImage,2)],'YData',[1 size(baseImage,1)]);
    
    image(cat(3,mat2gray(baseImage),mat2gray(xIn),zeros(size(baseImage))));
    axis off, axis image
    title('Aligned Image Overlay')
    
end
    
    


