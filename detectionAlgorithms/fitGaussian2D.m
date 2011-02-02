%FITGAUSSIAN2D Fit a 2-D Gaussian function to data in a square image window.
%    [prmVect prmStd C res J] = fitGaussian2D(data, prmVect, mode, options)
%
%    Symbols: xp : x-position
%             yp : y-position
%              A : amplitude
%              s : standard deviation
%              c : background
%
%    The origin is defined at the center of the input window.
%
%    Inputs:     data : 2-D image array
%             prmVect : parameter vector with order: [xp, yp, A, s, c]
%                mode : string that defines parameters to be optimized; any among 'xyasc'
%           {options} : vector [MaxIter TolFun TolX]; max. iterations, precision on f(x), precision on x
%
%    Outputs: prmVect : parameter vector
%              prmStd : parameter standard deviations
%                   C : covariance matrix
%                 res : residuals
%                   J : Jacobian
%
% Axis conventions: image processing, see meshgrid
% For Gaussian mixture fitting, use fitGaussianMixture2D()
%
% Example: [prmVect prmStd C res J] = fitGaussian2D(data, [0 0 max(data(:)) 1.5 min(data(:))], 'xyasc');

% (c) Francois Aguet & Sylvain Berlemont, 2011