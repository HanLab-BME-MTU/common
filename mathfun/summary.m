function s = summary(x, dim, stdflag)
% SUMMARY calculates the min, median, mean, max and standard deviation of a
% vector, matrix or N-D array.
% 
% SUMMARY(X) returns a structure, with the following fields:
% 
% MIN - For vectors, this is the minimum value of the elements in X.  For 
% matrices, this is a row vector containing the minimum value of each 
% column.  For N-D arrays, this is the minimum value of the elements along 
% the first non-singleton dimension of X.  See the function min for further
% information.
% 
% MEDIAN - For vectors, this is the median value of the elements in X.  For 
% matrices, this is a row vector containing the median value of each 
% column.  For N-D arrays, this is the median value of the elements along 
% the first non-singleton dimension of X.  See the function median for further
% information.
% 
% MEAN - For vectors, this is the mean value of the elements in X.  For 
% matrices, this is a row vector containing the mean value of each 
% column.  For N-D arrays, this is the mean value of the elements along 
% the first non-singleton dimension of X.  See the function mean for further
% information.
% 
% MAX - For vectors, this is the maximum value of the elements in X.  For 
% matrices, this is a row vector containing the maximum value of each 
% column.  For N-D arrays, this is the maximum value of the elements along 
% the first non-singleton dimension of X.  See the function max for further
% information.
% 
% STD - For vectors, this is the standard deviation of the elements in X.  For 
% matrices, this is a row vector containing the standard deviation of each 
% column.  For N-D arrays, this is the standard deviation of the elements along 
% the first non-singleton dimension of X.  See the function std for further
% information.
% 
% SIZE -  For vectors, matrices and N-D arrays, this is a vector of the
% sizes of each dimension of X. If X is a scalar, which MATLAB regards as a
% 1-by-1 array, SIZE(X) returns the vector [1 1]. 
% 
% 
% SUMMARY(X, DIM) calculates summary statistics along the dimension DIM of X.
% 
% 
% SUMMARY(X, DIM, STDFLAG) calculates summary statistics using a
% normalisation factor for standard deviation determined by STDFLAG. If
% STDFLAG is 0 (default), then the normalisation factor is n-1; if STDFLAG
% is 1, then the normalisation factor is n.  See the function std for
% further information.
% 
% 
% NOTE: Class type checking and error handling are conducted within the
% individual summary statistic calculation functions.
%
% 
% EXAMPLE: If X = [1 2 3; 2 4 6], then SUMMARY(X) returns
%            min: [1 2 3]
%         median: [1.5000 3 4.5000]
%           mean: [1.5000 3 4.5000]
%            max: [2 4 6]
%            std: [0.7071 1.4142 2.1213]
% 
% 
%   Class support for input X:
%      float: double, single
% 
% 
%   See also MIN, MEDIAN, MEAN, MAX, STD, SIZE.
% 
% 
% $ Author: Richie Cotton $     $ Date: 2006/03/16 $


% If the dimension was not specified, find the first nonsingleton
% dimension.
if nargin == 1 || dim == 0
    dim = find(size(x)~=1,1); 
    if isempty(dim)
        dim = 1;
    end
end

% If no flag was specified for standard deviation, default to 0.
if nargin < 3
   stdflag = 0; 
end

% Calculate summary statistics.
s.min = min(x, [], dim);
s.median = median(x, dim);
s.mean = mean(x, dim);
s.max = max(x, [], dim);
s.std = std(x, stdflag, dim);
s.size = size(x);