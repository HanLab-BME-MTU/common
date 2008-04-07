function [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data,dim,k)
%ROBUSTMEAN calculates mean and standard deviation discarding outliers
%
% SYNOPSIS [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data)
%
% INPUT    data : input data
%          dim  : (opt) dimension along which the mean is taken {1}
%          k    : (opt) #of sigmas at which to place cut-off {3}
%
% OUTPUT   finalMean : robust mean
%          stdSample : std of the data (divide by sqrt(n) to get std of the
%                      mean)
%          inlierIdx : index into data with the inliers 
%          outlierIdx: index into data with the outliers
%
% REMARKS  The code is based on (linear)LeastMedianSquares. It could be changed to
%          include weights 
%
% c: jonas, 04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% test input

if isempty(data) | ~all(isfinite(data)) %#ok<OR2>
    error('Please supply non-empty, finite data to robustMean')
end
if nargin<2 || isempty(dim)
    dim = 1;
end
if nargin < 3 || isempty(k)
    k = 3;
end


%========================
% LEAST MEDIAN SQUARES
%========================

% define magic numbers:
%k=3; %cut-off is roughly at 3 sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

% remember data size and reduced dataSize
dataSize = size(data);
reducedDataSize = dataSize;
reducedDataSize(dim) = 1;
% need this for later repmats
blowUpDataSize = dataSize./reducedDataSize;
% count how many relevant dimensions we have besides dim
realDimensions = length(find(dataSize>1));

% calc median - reduce dimension dim to length 1
medianData = median(data,dim);

% calculate statistics
res2 = (data-repmat(medianData,blowUpDataSize)).^2;
medRes2 = max(median(res2,dim),eps);

%testvalue to calculate weights
testValue=res2./repmat(magicNumber2*medRes2,blowUpDataSize);

if realDimensions == 1;
    %goodRows: weight 1, badRows: weight 0
    inlierIdx=find(testValue<=k^2);
    outlierIdx = find(testValue>k^2);

    % calculate std of the sample;
    if nargout > 1
        stdSample=sqrt(sum(res2(inlierIdx))/(length(inlierIdx)-4));
    end

    %====END LMS=========

    %======
    % MEAN
    %======

    finalMean = mean(data(inlierIdx));

else
    
    %goodRows: weight 1, badRows: weight 0
    inlierIdx=find(testValue<=k^2);
    outlierIdx=find(testValue > k^2);
    
    % mask outliers
    res2(outlierIdx) = NaN;
    % count inliers
    nInliers = sum(~isnan(res2),dim);

    % calculate std of the sample;
    if nargout > 1
        stdSample=sqrt(nansum(res2,dim)./(nInliers-4));
    end

    %====END LMS=========

    %======
    % MEAN
    %======
    
    data(outlierIdx) = NaN;
    finalMean = nanmean(data,dim);
end
