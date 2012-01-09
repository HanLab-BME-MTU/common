%[res, theta, nms, filterBank] = steerableDetector(img, M, sigma) performs edge/ridge detection through a generalization of Canny's alorithm based on steerable filters
%
% Inputs: 
%         img : input image
%           M : order of the filter, between 1 and 5
%             : Odd orders: edge detectors, M = 1 is equivalent to Canny's detector
%             : Even orders: ridge detectors
%               Higher orders provide better orientation selectivity and are less sensitive to noise,
%               at a small trade-off in computational cost.
%       sigma : standard deviation of the Gaussian kernel on which the filters are based
%
% Outputs: 
%         res : response to the filter
%       theta : orientation map
%         nms : non-maximum-suppressed response
%  filterBank : templates used in the computation of the steerable filter response
%
% For more information, see:
% M. Jacob et al., IEEE Trans. Image Proc., 26(8) 1007-1019, 2004.

% Francois Aguet, 07/12/2011 (last modified Nov 17, 2011).
% Adapted from SteerableJ package (2008).

function [res, theta, nms, filterBank] = steerableDetector(img, M, sigma) %#ok<STOUT,INUSD>