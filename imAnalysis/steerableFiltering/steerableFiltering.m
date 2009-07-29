function [R T NMS] = steerableFiltering(I, method, sigmaPSF) %#ok<STOUT,INUSD>
% Add some text.
%
% Input:
%
% I		input image.
%
% method	string designating which kind of steerable
%              filtering should be used. Choice are:
%              - '2ndGaussian' 2nd Derivative of a Gaussian
%              - 'UnserM1'     1st order Unser edge filter
%              - 'UnserM2'     2nd order Unser ridge filter
%
% sigma		standard deviation of the filter.
% 
% Output:
%
% R		filtering response.
%
% T		optimal angle.
%
% NMS		non-maximal suppression (optional).

end
