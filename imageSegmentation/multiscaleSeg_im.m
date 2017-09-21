function refinedMask = multiscaleSeg_im(im, varargin)
% multiscaleSeg Segment a single cell image by combining segmentation
% results obtained at multiple smoothing scales. Since it requires only one
% tuning parameters (tightness) and ‘tightness’=0.5 works well for many cases, 
% it achieves almost automatic segmentation.
%
% Input: an image
% Output: refined mask
%
% 2017/05/29, Jungsik Noh


ip = inputParser;
ip.addParameter('type', 'middle');
ip.addParameter('tightness', -1);
ip.addParameter('numModels', -1);

ip.parse(varargin{:});
p = ip.Results;


%% type

switch p.type
    case 'middle'
    case 'tight'
    case 'minmax'
    otherwise
        error('Type should be either of "middle", "tight" or "minmax".') 
end


%%  Parameters
% figFlag = 'off';

sigmas = [0 0.66 1 1.66 2.5 4];  % unit: pixel (common scales for xxx by xxx size confocal images)

p.MinimumSize = 100;        
p.ObjectNumber = 1;
p.FillHoles = 1;

%disp(p)

numModels = numel(sigmas)*2*2;

maskingResultArr = nan(size(im, 1), size(im, 2), numModels);


%%  Gaussian filtering, method, refinement

for k = 1:numel(sigmas)
    %disp(['GaussFilterSigma: ', num2str(sigmas(k))])

    currImage = im;

% GaussianFiltering
if sigmas(k) > 0
    currImage = filterGauss2D(currImage, sigmas(k));
end


% minmax
try
    currThresh = thresholdFluorescenceImage(currImage); 
    currMask1 = (currImage >= currThresh);
    
    p.ClosureRadius = 1;
    refinedMask1 = maskRefinementCoreFunc(currMask1, p);
    maskingResultArr(:,:, 4*(k-1)+1) = refinedMask1;

    p.ClosureRadius = 3;
    refinedMask1 = maskRefinementCoreFunc(currMask1, p);
    maskingResultArr(:,:, 4*(k-1)+2) = refinedMask1;
    
catch
    disp(['GaussFilterSigma: ', num2str(sigmas(k))])
    disp('Error in minmax thresholding')
    maskingResultArr(:,:, 4*(k-1)+1) = nan(size(currImage, 1), size(currImage, 2));
    maskingResultArr(:,:, 4*(k-1)+2) = nan(size(currImage, 1), size(currImage, 2));
end
    

% Rosin
currThresh = thresholdRosin(currImage); 
currMask1 = (currImage >= currThresh); 

p.ClosureRadius = 1;
refinedMask1 = maskRefinementCoreFunc(currMask1, p);
maskingResultArr(:,:, 4*(k-1)+3) = refinedMask1;

p.ClosureRadius = 3;
refinedMask1 = maskRefinementCoreFunc(currMask1, p);
maskingResultArr(:,:, 4*(k-1)+4) = refinedMask1;

end



%% sum maskings from multiple methods

res = sum(maskingResultArr(:,:,1:end), 3, 'omitnan');
tab = tabulate(res(:));
 tabulate(res(:)) 

%figure('Visible', figFlag);
%imagesc(res)
%colormap(jet);colorbar
%pause(1)

val = tab(:,1);
counts = tab(:,2);

[~, ind] = sort(counts, 'descend');

a = val(ind(1));
b = val(ind(2));
backgroundth = min(a, b);
maskth = max(a, b);


%% ensemble method

if ((p.tightness < 0) && (p.numModels < 0))

switch p.type
    case 'middle'
    % ensemble method: median number of models
%        halfEffNumModel = sum(~isnan(maskingResultArr(1,1,:)))/2;
        midNumModel0 = (backgroundth+1+maskth)/2;
        midNumModel = midNumModel0;
%        midNumModel = max(midNumModel0, halfEffNumModel);  
        % middle of two peaks should be greater than 50% voting. (xx)

        res0 = (res > midNumModel);
        mnum_smallest = maskth;
        mnum_biggest = backgroundth + 1;
        tightness_interp = interp1([mnum_biggest, mnum_smallest], [0, 1], midNumModel);        
        disp('Threshold of votes:'); 
        disp([num2str(midNumModel), ' (tightness: ', num2str(tightness_interp), ')'])

    case 'tight'
    % ensemble: tight        
        res0 = (res >= maskth);
        mnum_smallest = maskth;
        mnum_biggest = backgroundth + 1;
        tightness_interp = interp1([mnum_biggest, mnum_smallest], [0, 1], maskth);        
        disp('Threshold of votes:'); 
        disp([num2str(maskth), ' (tightness: ', num2str(tightness_interp), ')'])
        
    case 'minmax'
        mnumint = max(ind(1), ind(2))-1:-1:min(ind(1), ind(2))+1;
        tmp = counts(mnumint);
        tmp2 = tmp(2:end) + tmp(1:end-1);
        dif0 = reshape( diff( tmp2 ), [], 1);
        delta = find(([dif0; 1] > 0), 1);
        localMinIndex0 = max(ind(1), ind(2)) - delta - 1;
        localMinIndex = max(localMinIndex0, 1);
        minmaxNumModel = max(val(localMinIndex), backgroundth+1);
        
        res0 = (res > minmaxNumModel);
        mnum_smallest = maskth;
        mnum_biggest = backgroundth + 1;
        tightness_interp = interp1([mnum_biggest, mnum_smallest], [0, 1], minmaxNumModel);        
        disp('Threshold of votes:'); 
        disp([num2str(minmaxNumModel), ' (tightness: ', num2str(tightness_interp), ')'])
end

elseif p.tightness >= 0
    if p.tightness > 1
        error('Tightness should range from 0 to 1.')
    end
    
    mnum_smallest = maskth;
    mnum_biggest = backgroundth + 1;
    mnum_interp = interp1([0, 1], [mnum_biggest, mnum_smallest], p.tightness);
    tightnessNumModel = round(mnum_interp);
    
    res0 = (res > tightnessNumModel);
    disp('Threshold of votes:'); 
    disp([num2str(tightnessNumModel), ' (tightness: ', num2str(p.tightness), ')'])

else
    res0 = (res >= p.numModels);
    disp('Threshold of votes:'); 
    disp(num2str(p.numModels))
    
end

 

%% final refinement

    p.ClosureRadius = 1;
    refinedMask = maskRefinementCoreFunc(res0, p);
 
%bdd = bwboundaries(refinedMask);
%bdd1 = bdd{1};
%figure('Visible', figFlag);
%imshow(im, [])
%hold on
%plot(bdd1(:,2), bdd1(:,1), 'r');
%hold off

%tmp = (res0~=refinedMask);
%disp(sum(tmp(:)))


%disp('for i=1:30; close(figure(i)); end')

end

