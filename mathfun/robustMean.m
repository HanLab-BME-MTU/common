function [finalMean, stdSample, inlierIdx] = robustMean(data)
%ROBUSTMEAN calculates mean and standard deviation discarding outliers
%
% SYNOPSIS [finalMean, stdSample, inlierIdx] = robustMean(data)
%
% INPUT    data : input data
%          
% OUTPUT   finalMean : robust mean
%          stdSample : std of the data (divide by sqrt(n) to get std of the
%                      mean)
%          inlierIdx : index into data with the inliers (recover the others 
%                      with findMissingIndices)
%
% REMARKS  The code is based on (linear)LeastMedianSquares. It could be changed to
%          include weights and/or application of the mean along other
%          dimensions
%
% c: jonas, 04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% test input

if isempty(data) | ~all(isfinite(data))
    error('Please supply non-empty, finite data to robustMean')
end


%========================
% LEAST MEDIAN SQUARES
%========================

% define magic numbers:
k=3; %value important for calculation of sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

% prepare data
data = data(:);

% calc median
medianData = median(data);

% calculate statistics
res2 = (data-medianData).^2;
medRes2 = median(res2);

%testvalue to calculate weights
testValue=res2/(magicNumber2*medRes2);

%goodRows: weight 1, badRows: weight 0
inlierIdx=find(testValue<=k^2);

% calculate std of the sample;
if nargout > 1
    stdSample=sqrt(sum(res2(inlierIdx))/(length(inlierIdx)-4));
end

%====END LMS=========

%======
% MEAN
%======

finalMean = mean(data(inlierIdx));
