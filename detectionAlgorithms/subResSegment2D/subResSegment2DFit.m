function [F J] = subResSegment2DFit(varSegmentParams, I, fixSegmentParams, paramSelector)
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
%                   and pVar is the number of parameters being optimized.
%                   This is a sub-set of the whole set of segment
%                   parameters (i.e. xC, yX, amp, sigma, l, theta). pVar is
%                   subject to the following requirements:
%                   - pVar == nnz(paramSelector)
%                   - pVar + size(fixSegmentParams,2) == numel(paramSelector)
%
% I                  image
%
% fixSegmetnParams  (optional) a n by pFixed matrix where n is the number
%                   of segments and pFix is the number of parameters not
%                   being optimized. It is a sub-set of the whole set of
%                   segment parameters (i.e. xC, yX, amp, sigma, l, theta).
%                   pFix is subject to the following requirements:
%                   - pFix == ~nnz(paramSelector)
%                   - size(varSegmentParams,2) + pFix == numel(paramSelector)
%
% paramSelector     (optional) a logical vector of size 6 where
%                   paramSelector(i) == true means that ith segment
%                   parameter will be optimized. By default all parameters
%                   are involved in the optimization. The number of true
%                   elements must match the number of columns of
%                   varSegmentParams ; the number of false elements must
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

if nargin ~= 2 || nargin ~= 4
    error('Invalid number of input arguments.');
end

[n pVar] = size(varSegmentParams);

if pVar > p
    error('Too many segment parameters.');
end

if nargin > 2
    if ~islogical(paramSelector)
        error('paramSelector must be logical.');
    end
    
    if numel(paramSelector) ~= p
        error('paramSelector must be of size 6.');
    end
    
    pFix = size(fixSegmentParams,2);
    
    if pVar + pFix ~= numel(paramSelector) || pVar ~= nnz(paramSelector)
        error('1st and 3rd parameter size do not match paramSelector.');
    end
    
    if size(fixSegmentParams,1) ~= n
        error('1st and 3rd parameter size differ.');
    end
else
    fixSegmentParams = [];
    paramSelector = true(p,1);
end

% Merge varSegmentParams and fixSegmentParams
segmentParams = zeros(n,p);
paramSet = 1:p;
segmentParams(:,paramSet(paramSelector == true)) = varSegmentParams;
segmentParams(:,paramSet(paramSelector == false)) = fixSegmentParams;

[nrows ncols] = size(I);

m = nrows * ncols;
    
xRange = cell(n, 1);
yRange = cell(n, 1);
nzIdx = cell(n,1);

F = zeros(size(I));

for i = 1:n
    xC = segmentParams(i,1);
    yC = segmentParams(i,2);
    amp = segmentParams(i,3);
    sigma = segmentParams(i,4);
    l = segmentParams(i,5);
    theta = segmentParams(i,6);
    
    [xRange{i},yRange{i},nzIdx{i}] = ...
        subResSegment2DSupport(xC,yC,sigma,l,theta,[nrows,ncols]);

    S = subResSegment2D(xRange{i}-xC,yRange{i}-yC,amp,sigma,l,theta,nzIdx{i});

    F(yRange{i},xRange{i}) = F(yRange{i},xRange{i}) + S;
end

F = reshape(F - I, m, 1);

if nargout > 1
    indPixels = cell(n,1);
    indParams = cell(n,1);
    val = cell(n,1);    
    
    for i = 1:n
        xC = segmentParams(i,1);
        yC = segmentParams(i,2);
        amp = segmentParams(i,3);
        sigma = segmentParams(i,4);
        l = segmentParams(i,5);
        theta = segmentParams(i,6);
        
        % Compute all partial derivatives of F against segment parameters
        dF = subResSegment2DJacobian(xRange{i}-xC,yRange{i}-yC,amp,sigma,...
            l,theta,nzIdx{i},paramSelector);
        
        [X Y] = meshgrid(xRange{i}, yRange{i});
        ind = sub2ind([nrows ncols], Y(:), X(:));
        
        indPixels{i} = repmat(ind, p, 1);
        indParams{i} = cell2mat(arrayfun(@(k) ones(numel(ind), 1) * i +...
            k * n, (0:p-1)', 'UniformOutput', false));

        val{i} = dF(:);
    end
    
    indPixels = vertcat(indPixels{:});
    indParams = vertcat(indParams{:});
    val = vertcat(val{:});
    J = sparse(indPixels, indParams, val, m, n * pVar, length(val));
end