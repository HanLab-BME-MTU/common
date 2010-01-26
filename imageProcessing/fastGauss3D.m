function out=fastGauss3D(img,sigma,fSze,correctBorder,filterMask,reduceNanImage)
% fastGauss3D	apply a 2/3 dimensional gauss filter
%
%    SYNOPSIS out=fastGauss3D(img,sigma,fSze,correctBorder,filterMask)
%
%    INPUT:   img      2 or 3-dimensional data
%             sigma    of gauss filter. Supply single sigma or sigma for
%                      every dimension
%             fSze     (optional) size of the gauss mask [sizeX sizeY sizeZ]
%                          (odd size required for symmetric mask!)
%                       If empty, odd mask with +/- 4*sigma is used.
%             correctBorder (optional) if 1, border effects of the filtering are
%                              lessened. If 2, old version of correction is
%                              used. Default: 1. Note: correctBorder = 2
%                              requires 3D image.
%             filterMask    (optional) supply your own mask. In this case,
%                              sigma and fSze don't have to be supplied. If
%                              the filter mask is a cell, the data will be
%                              sequentially filtered by the contents of the
%                              cell. In that case, fSze is required.
%             reduceNanImage (optional) if 1, nan-containing areas are
%                              clipped before filtering to increase speed.
%                              Default: 0
%
%    OUTPUT:  out      filtered data
%
% c: 13/03/01 dT
% revamped by jonas

%===============
% check input
%===============
if nargin < 1 || isempty(img);
    error('Please pass nonempty image to fastGauss3D')
end
% set default correctBorder
if nargin < 4 || isempty(correctBorder)
    correctBorder = 1;
end

% read dimensionality
dims = sum(size(img)>0);

% check for filterMask
if nargin < 5 || isempty(filterMask)
    % in this case, we need fSze etc
    if nargin < 2 || isempty(sigma) || ~any(length(sigma)==[1 dims])
        error('please supply nonempty sigma of correct dimensionality')
    end
    if length(sigma) == 1
        sigma = sigma * ones(1,dims);
    end
    % check for filterSize
    if nargin < 3 || isempty(fSze)
        fSze = roundOddOrEven(sigma(1:dims)*4,'odd','inf');
    end
    % create gaussMask
    switch dims
        case 2
            filterMask=GaussMask2D(sigma,fSze,[],1,[],[],1);
        case 3
            filterMask=GaussMask3D(sigma,fSze,[],1,[],[],1);
    end
else
    if isempty(fSze)
        if iscell(filterMask)
            error('if you supply a separated filter, you will need to supply the filter size!')
        else
            fSze = size(filterMask);
        end
    end
end

if nargin < 6 || isempty(reduceNanImage)
    reduceNanImage = false;
end


% add border to image
convnOpt = 'same';
nanMask = []; % mask indicating all the NaNs in the original image;
if correctBorder == 2
    if dims < 3
        warning('FASTGAUSS3D:WRONGDIMENSION',...
            'Cannot use old correctBorder if not 3D image. Using new correctBorder instead')
        correctBorder = 1;
    else
        addBorderOld;
        
        halFsze = zeros(1,3);
        halFsze(1:length(fSze)) = floor(fSze/2); %half filter size
        
        
        img = padarrayXT(img, [halFsze halFsze halFsze], 'symmetric');
        convnOpt = 'valid';
    end
end
if correctBorder == 1
    %correct for border effects. Create nanMask first. It's a logical mask
    %that takes up 1/8th the space of the original. As long as more than
    %1/8th of the image contains Nan, it's smaller than a list of indices
    nanMask = isnan(img);
    fullMask = [];
    if reduceNanImage && any(nanMask(:))
        goodRCZcell = cell(dims,1);
        nan2 = all(nanMask,3);
        goodRCZcell{1} = ~all(nan2,2);
        goodRCZcell{2} = ~all(nan2,1);
        clear nan2
        if dims > 2
            goodRCZcell{3} = ~all(all(nanMask,1),2);
        end
        nanRatio = 1-prod(cellfun(@(x)(sum(x)),goodRCZcell))/numel(nanMask);
        if nanRatio > 0.05 %-- from a little bit of testing it looks like you get about nanRatio*0.4 reduction in time
            img = img(goodRCZcell{:});
            fullMask = nanMask;
            nanMask = nanMask(goodRCZcell{:});
        end
    end
    
    % pass nanMask so that we don't need to recalc again
    img = padarrayXT(img, floor(fSze/2)*ones(1,dims), 'symmetric');
    convnOpt = 'valid';
end

switch dims
    case 2
        % Convolve matrices
        if iscell(filterMask)
            for i=1:length(filterMask)
                if length(filterMask{i}) > 1
                    img = conv2(img,filterMask{i},convnOpt);
                end
            end
            out = img;
        else
            out=conv2(img,filterMask,convnOpt);
        end
    case 3
        % Convolve matrices
        if iscell(filterMask)
            for i=1:length(filterMask)
                if length(filterMask{i}) > 1
                    img = convn(img,filterMask{i},convnOpt);
                end
            end
            out = img;
        else
            out=convn(img,filterMask,convnOpt);
        end
end


if ~isempty(nanMask) && reduceNanImage && ~isempty(fullMask)
    % rebuild full image
    img = out;
    out = NaN(size(fullMask));
    out(~fullMask) = img(~nanMask);
else
    % simply ensure that NaNs are in the right place
    out(nanMask) = NaN;
end