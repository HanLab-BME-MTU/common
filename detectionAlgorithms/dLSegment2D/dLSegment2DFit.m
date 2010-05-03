function [F J] = subResSegment2DFit(x, I, sigmaPSF)
% This function computes I - Im, where I is an image and Im is an image
% model defined by the sum of sub-resolution 2D segment model caracterize
% by params x (see subResSegment2DImageModel() function for more details).
% This function intends to be used as a function handler in lsqnonlin,
% lsqlin, etc. optimization functions.
%
% [F J] = subResSegment2DFit(x, I, sigmaPSF)
%
% Parameters:
% x            a nx6 matrix where n is the number of segments and their
%              parameters, i.e. xC, yC, A, l, t are stored column-wise
%              (see subResSegment2D.m for more details)
%
% I            image
%
% sigmaPSF     half width of the gaussian PSF model
%
% Output:
% F            F = model - I
%
% J            J is an MxN matrix where M = numel(I) and N = numel(x). This
%              corresponds to the jacobian of I against x (see lsqnonlin for
%              mode details)

[n p] = size(x);

if p ~= 5
    error('Invalid number of segment parameters.');
end

[nrows ncols] = size(I);

m = nrows * ncols;
    
xRange = cell(n, 1);
yRange = cell(n, 1);
nzIdx = cell(n,1);

F = zeros(size(I));

for i = 1:n
    xC = x(i,1);
    yC = x(i,2);
    A = x(i,3);
    l = x(i,4);
    t = x(i,5);
    
    [xRange{i},yRange{i},nzIdx{i}] = subResSegment2DSupport(xC,yC,sigmaPSF,...
        l,t,[nrows,ncols]);

    S = subResSegment2D(xRange{i}-xC,yRange{i}-yC,A,sigmaPSF,l,t,nzIdx{i});

    F(yRange{i},xRange{i}) = F(yRange{i},xRange{i}) + S;
end

F = reshape(F - I, m, 1);

if nargout > 1
    indPixels = cell(n,1);
    indParams = cell(n,1);
    val = cell(n,1);    
    
    for i = 1:n
        xC = x(i,1);
        yC = x(i,2);
        A = x(i,3);
        l = x(i,4);
        t = x(i,5);
        
        % Compute all partial derivatives of F against segment parameters
        [dFdXc, dFdYc, dFdA, dFds, dFdl, dFdt] = ...
            subResSegment2DJacobian(xRange{i}-xC, yRange{i}-yC, A, sigmaPSF, ...
            l, t,nzIdx{i}); %#ok<ASGLU>
        
        [X Y] = meshgrid(xRange{i}, yRange{i});
        ind = sub2ind([nrows ncols], Y(:), X(:));
        
        indPixels{i} = repmat(ind, p, 1);
        indParams{i} = cell2mat(arrayfun(@(k) ones(numel(ind), 1) * i + ...
            k * n, (0:p-1)', 'UniformOutput', false));
        val{i} = vertcat(dFdXc(:), dFdYc(:), dFdA(:), dFdl(:), dFdt(:));
    end
    
    indPixels = vertcat(indPixels{:});
    indParams = vertcat(indParams{:});
    val = vertcat(val{:});
    J = sparse(indPixels, indParams, val, m, n * p, length(val));
end