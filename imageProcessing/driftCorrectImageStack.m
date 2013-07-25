function imStackCorr = driftCorrectImageStack(imStack,varargin)
%DRIFTCORRECTIMAGESTACK performs drift correction / image stabilization / intra-stack registration on the input image stack
%
% [imStackCorr,tForms] = driftCorrectImageStack(imStack)
%
% This function attempts to correct for drift/misalignment in the input
% image stack in a two-step approach. First the images are aligned by
% maximum cross-correlation (pixel-level translation only), and then this
% intiial alignment is (optionally) refined using imregister (sub-pixel,
% supports translation in addition to rigid, affine, etc.)
%
% Input:
%
%   imStack - 3D stack of 2D images. (Though it wouldn't be too much work
%   to generalize this to 4D stacks of 3D... feel free to do so!!)
%
%   Also, a bunch of parameters that you'll have to open the function to
%   see descriptions of because I'm too lazy to re-document here.
%
%Hunter Elliott
%7/2013


%% -------------- Input ------------ %%

ip = inputParser;
ip.addRequired('imStack')
ip.addParamValue('ReferenceImage',-1,@(x)(isscalar(x) && isequal(x,round(x)) && x ~= 0 ));%Image to correct to. Negative numbers are relative e.g. -1 = previos image
ip.addParamValue('DoRefinement',true,@islogical);%Whether to run imregister after the correlation-based translation to improve the registration
ip.addParamValue('TransformType','translation');%See imregister.m
ip.addParamValue('Modality','monomodal')%See imregister.m
ip.addParamValue('FillValue',0,@isscalar);%Value to fill in pixels which are empty due to correction
ip.parse(imStack,varargin{:});

p = ip.Results;

%% ----------- Init --------- %%

ogClass = class(imStack);
imSize = size(imStack);

imStackCorr = imStack;

%Get indices of reference images and corresponding images to correct
if p.ReferenceImage > 0
    %If a constant reference image was specified
    refInd = ones(imSize(end),1) * p.ReferenceImage;
    corrInd = [1:p.ReferenceImage-1 p.ReferenceImage+1:imSize(end)];
else
    refInd = 1:imSize(end)+p.ReferenceImage;
    corrInd = -p.ReferenceImage+1:imSize(end);    
end

nCorr = numel(corrInd);

[optimizer,metric] = imregconfig(p.Modality);

showPlots = true;

%% ----------- Correction --------- %%

for j = 1:nCorr
       
    currRef = double(imStackCorr(:,:,refInd(j)));    
        
    %Get cross-orrelation between images to do an initial translation since imregister is slow (optimization
    %based)    
    imXcorr = convnfft(currRef,rot90(imStack(:,:,corrInd(j)),2));
    %Find max correlation and shift accordingly
    [~,iMax] = max(imXcorr(:));    
    [maxCorr(1), maxCorr(2)] = ind2sub(size(imXcorr),iMax);        
    [Xi,Yi] = meshgrid((1:imSize(2))+imSize(2)-maxCorr(2),(1:imSize(1))+imSize(1)-maxCorr(1));
    currCorr = interp2(double(imStack(:,:,corrInd(j))),Xi,Yi,'nearest');%Lazy way to do shifting...
    currCorr(isnan(currCorr)) = p.FillValue;
    
    if showPlots
        fsFigure(.75);
        subplot(1,2,1)
        imshow(cat(3,mat2gray(currRef),mat2gray(imStack(:,:,corrInd(j))),zeros(imSize(1:2))))        
        subplot(1,2,2)
        imshow(cat(3,mat2gray(currRef),mat2gray(currCorr),zeros(imSize(1:2))))        
    end
    
    %Optionally improve the registration
    if p.DoRefinement
        imStackCorr(:,:,corrInd(j)) = imregister(currCorr,currRef,p.TransformType,optimizer,metric);
    else        
        imStackCorr(:,:,corrInd(j)) = currCorr;
    end
    
    if showPlots
        fsFigure(.75);
        subplot(1,2,1)
        imshow(cat(3,mat2gray(currRef),mat2gray(currCorr),zeros(imSize(1:2))))        
        subplot(1,2,2)
        imshow(cat(3,mat2gray(currRef),mat2gray(imStackCorr(:,:,corrInd(j))),zeros(imSize(1:2))))
    end
    
    
end


