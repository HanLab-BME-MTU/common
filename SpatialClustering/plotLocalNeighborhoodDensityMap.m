function [C, H] = plotLocalNeighborhoodDensityMap(mpm,imsize,varargin)

ip = inputParser;
ip.CaseSensitive = false;
%simulation params
ip.addRequired('mpm', @isnumeric);
ip.addRequired('imsize', @isnumeric);
ip.addParamValue('dist', 1:20, @isnumeric);
ip.addParamValue('scale', 1, @isscalar);
ip.addParamValue('lr', [], @isnumeric);
ip.parse(mpm, imsize, varargin{:});
dist = ip.Results.dist;
lr = ip.Results.lr;

%MAKE MASK
imsizS = [imsize(2) imsize(1)];
%[areamask] = makeCellMaskDetections(mpm,closureRadius,...
%dilationRadius,doFill,imsize,plotMask,[]);
%CALCULATE NORMALIZED AREA FROM MASK
%normArea = bwarea(areamask);
% CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
%corrFacMat  = makeCorrFactorMatrix(imsizS, dist, 10, areamask');


%Calculate L-ripley
if isempty(lr)
    lr = nan(length(dist),size(mpm,1));
    for i = 1:length(mpm)
        [~,lr(:,i)]=RipleysKfunction(mpm([1:i-1 i+1:length(mpm)],:),mpm(i,:),imsizS,dist); %,corrFacMat,normArea);
    end
end

%make map out of L values for each point at a given scale (i.e., for each
%point assign it's L(d) for a given d)
M = zeros(imsize);
for i = 1:size(lr,2)
    M(max(1,round(mpm(i,2))),max(1,round(mpm(i,1)))) = lr(ip.Results.scale,i);
end


%add random points to make sure plot fills out nicely
percentPoints = 0.5;
mpmRandSample = repmat(imsize,0.25*imsize(1)*imsize(2),1).*...
    rand(0.25*imsize(1)*imsize(2),2);
lrRandSample = nan(length(dist),size(mpm,1));
for i = 1:length(mpmRandSample)
    [~,lrRandSample(:,i)]=RipleysKfunction(mpm,mpmRandSample(i,:),imsizS,dist); %,corrFacMat,normArea);
end
for i = 1:size(lrRandSample,2)
    M(max(1,round(mpmRandSample(i,2))),max(1,round(mpmRandSample(i,1)))) = lrRandSample(ip.Results.scale,i);
end


%plot contours
V = linspace(min(M(:)),max(M(:)),10);
figure
[C H] = contourf(M);
clabel(C,H);
hold on
plot(mpm(:,1),mpm(:,2),'m.')