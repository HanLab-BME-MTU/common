function [R T NMS FBANK] = steerableFiltering(I, M, sigmaPSF) %#ok<STOUT,INUSD>

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used. Use steerableDetector instead.\n']);

[R T NMS FBANK] = steerableDetector (I, M, sigmaPSF);
% [R T NMS FBANK] = steerableFiltering(I, M, sigmaPSF)
%
% This function implements the article:
% Jacob M, Unser M. "Design of steerable filters for feature detection
% using canny-like criteria" IEEE Trans Pattern Anal Mach Intell. 2004
% Aug;26(8):1007-19.
%
% Input:
%
% I		    input image.
%
% M	        filter order (M = 1, 2, 3 or 4). Use an even order for ridge
%           detection and an odd order for edge detection (M = 1
%           corresponds to Canny edge detector).
%
% sigma		standard deviation of the filter.
% 
% Output:
%
% R		    filtering response.
%
% T		    optimal angle.
%
% NMS		non-maximal suppression.
%
% FBANK     the image convolved by every basis filter from which the Mth
%           order filter is made of. The number of basis filter is defined
%           by the recursive formula
%           nfilter = acc(M),
%           where acc(i) = (i*(i+3)/2) - acc(i-1)
%                 acc(0) = 0

end
