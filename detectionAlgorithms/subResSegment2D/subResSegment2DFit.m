function [F J] = subResSegment2DFit(varSegmentParams, I, avgBkg, fixSegmentParams, paramSelector)
% This function computes I - Im, where I is an image and Im is an image
% model defined by the sum of sub-resolution 2D segment model caracterize
% by params (the first and two last parameters). This function intends to
% be used as a function handler in lsqnonlin, lsqlin, etc. optimization
% functions. The initialization of segmentParams can be obtained by
% subResSegment2DInit function.
%
% [F J] = subResSegment2DFit(varSegmentParams, I, fixSegmentParams, paramSelector)
%
% Parameters:
% varSegmentParams  a n by pVar matrix where n is the number of segments
%                   and pVar is the number of parameters being optimized
%                   (excepted background value treated as a special case).
%                   This variable is a sub-set of the whole set of segment
%                   parameters (i.e. x, y, amp, sigma, l, theta). pVar is
%                   subject to the following requirements:
%                   - pVar == nnz(paramSelector) - paramSelector(7)
%                   - pVar + size(fixSegmentParams,2) == 6
%
% I                  image
%
% fixSegmetnParams  (optional) a n by pFix matrix where n is the number
%                   of segments and pFix is the number of parameters not
%                   being optimized. It is a sub-set of the whole set of
%                   segment parameters (i.e. xC, yX, amp, sigma, l, theta).
%                   pFix is subject to the following requirements:
%                   - pFix == ~nnz(paramSelector)
%                   - size(varSegmentParams,2) + pFix == numel(paramSelector)
%
% paramSelector     (optional) a logical vector of size 7 where
%                   paramSelector(i) == true means that ith segment
%                   parameter will be optimized. By default all parameters
%                   are involved in the optimization. The number of true
%                   elements must be equal to the number of columns of
%                   varSegmentParams + 1 ; the number of false elements must
%                   match the number of columnds of fixSegmentParams.
%                   
% Output:
% F            F = model - I
%
% J            J is an MxN sparse matrix where M = numel(I) and N =
%              numel(varSegmentParams). This corresponds to the jacobian of
%              I against segmentParams (see lsqnonlin for mode details).
%
% Sylvain Berlemont, 2010

p = 6;

if nargin ~= 3 && nargin ~= 5
    error('Invalid number of input arguments.');
end

[n pVar] = size(varSegmentParams);

if pVar > p
    error('Too many segment parameters.');
end

if nargin > 3
    if ~islogical(paramSelector)
        error('paramSelector must be logical.');
    end
    
    if numel(paramSelector) ~= p
        error(['paramSelector must be of size ' num2str(p)]);
    end
    
    pFix = size(fixSegmentParams,2);
    
    if pVar + pFix ~= numel(paramSelector) || pVar ~= nnz(paramSelector)
        error('1st and 3rd parameter size do not match paramSelector.');
    end
    
    if size(fixSegmentParams,1) ~= n
        error('1st and 3rd parameter size differ.');
    end
else
    fixSegmentParams = zeros(n,0);
    paramSelector = true(p,1);
end

% Merge varSegmentParams and fixSegmentParams
segmentParams = zeros(n,p);
paramSet = 1:p;
segmentParams(:,paramSet(paramSelector)) = varSegmentParams;
segmentParams(:,paramSet(~paramSelector)) = fixSegmentParams;

[nrows ncols] = size(I);

m = nrows * ncols;

% Get the image model
[Im xRange yRange nzIdx] = subResSegment2DImageModel(segmentParams, avgBkg, [nrows ncols]);

% Compute F
F = reshape(Im - I, m, 1);

if nargout > 1
    indPixels = cell(n,1);
    indParams = cell(n,1);
    val = cell(n,1);    
    
    for i = 1:n
        x = segmentParams(i,1);
        y = segmentParams(i,2);
        amp = segmentParams(i,3);
        sigma = segmentParams(i,4);
        l = segmentParams(i,5);
        theta = segmentParams(i,6);
        
        % Compute all partial derivatives of F against segment parameters
        dF = subResSegment2DJacobian(x,y,amp,sigma,...
            l,theta,xRange{i},yRange{i},nzIdx{i},paramSelector);
        
        [X Y] = meshgrid(xRange{i}, yRange{i});
        ind = sub2ind([nrows ncols], Y(:), X(:));
        
        indPixels{i} = repmat(ind, pVar, 1);
        indParams{i} = cell2mat(arrayfun(@(k) ones(numel(ind), 1) * i +...
            k * n, (0:pVar-1)', 'UniformOutput', false));

        val{i} = dF(:);
    end
    
    indPixels = vertcat(indPixels{:});
    indParams = vertcat(indParams{:});
    val = vertcat(val{:});
    J = sparse(indPixels, indParams, val, m, n * pVar, length(val));
end