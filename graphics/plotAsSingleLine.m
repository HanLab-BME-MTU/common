function [ h ] = plotAsSingleLine( varargin )
%plotAsSingleLine is a drop-in replacement for plot, but tries to create as
%few graphics objects by only creating one line object per LineSpec.
%
% The input / output scheme should be exactly like plot, except that the
% number handles output will correspond to the number of LineSpecs
%
% See also plot

% Mark Kittisopikul, January 2015

%% 1. Check if the first argument is an axis handle.
% If so store and shift inputs.

inIdx = 1;
if(isscalar(varargin{inIdx}) && ishandle(varargin{inIdx}))
    ax = varargin{inIdx};
    inIdx = inIdx + 1;
    varargin = varargin(2:end);
else
    ax = gca;
end

%% 2. Collect information about input type and dimensions

isXYdata = cellfun(@isnumeric,varargin);
nRows = cellfun('size',varargin,1);
nCols = cellfun('size',varargin,2);
numEl = nRows.*nCols;
% isLineSpec = cellfun(@ischar,varargin);
isVector = nRows == 1 | nCols == 1;
isRowVector = nRows == 1;

%% 3. If we are only given Y data, then plot that and return

if(sum(isXYdata) == 1)
    % only Y data is given
    [X,Y] = joinColumns(NaN, repmat( (1:nRows(1))' ,1,nCols(1)) ,varargin{inIdx});
    h = plot(X,Y,varargin{inIdx+1:end});
    return;
end

%% 4. Otherwise, X and Y data must come in pairs.
% Determine which input is X and which is Y
isXYpairX = [isXYdata(1:end-1) & isXYdata(2:end) false];
isXYpairX(isXYdata) = mod(1:sum(isXYdata),2);
isXYpairY = [false isXYpairX(1:end-1)];
isXYpair = isXYpairX | isXYpairY;

%% 5. Turn all vectors into column vectors.
isRowVector = isRowVector & isXYpair;
varargin(isRowVector) = cellfun(@transpose,varargin(isRowVector),'UniformOutput',false);

%% 6. Transpose matrices if the vector and matrix dimensions do not agree

%if the number of rows in the column vector does not match the number of
%rows in it's partner, then transpose the partner matrix

%6a: first if X is a vector, then tranpose the Ys
filterX = isVector & isXYpairX;
filterY = [false filterX(1:end-1)];
filter = filterY;
filter(filterY) = numEl(filterX) ~= nRows(filterY);
varargin(filter) = cellfun(@transpose,varargin(filter),'UniformOutput',false);

%6b: next if Y is a vector, then transpose the Xs
filterY = isVector & isXYpairY;
filterX = [filterY(2:end) false];
filter = filterX;
filter(filterX) = nRows(filterX) ~= numEl(filterY);
varargin(filter) = cellfun(@transpose,varargin(filter),'UniformOutput',false);

%% 7. Expand vectors to match the size of matrices
% Try to make X and Y the same size
[varargin(isXYpairX),varargin(isXYpairY)] = ...
    cellfun(@(x,y) deal( ...
    repmat(x,ceil(size(y)./size(x))), ...
    repmat(y,ceil(size(x)./size(y)))), ...
    varargin(isXYpairX), ...
    varargin(isXYpairY), ...
    'UniformOutput',false);
% no guarantee they are the same size, but plot will catch dimension
% mismatch
% numEl = cellfun('prodofsize',varargin);
% assert(all(numEl(isXYpairX) == numEl(isXYpairY)),'Vectors must be the same length');

%% 8. Pack lines into a single column per linespec
[varargin{isXYpair}] = joinColumns(NaN,varargin{isXYpair});

%% 9. Forward data to plot
h = plot(ax,varargin{:});


end

