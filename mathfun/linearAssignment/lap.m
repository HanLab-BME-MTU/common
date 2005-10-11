function [x, y] = lap(cc, NONLINK_MARKER, extendedTesting, augmentCC)
%LAP solves the linear assignment problem for a given cost matrix
%
% A linear assignment tries to establish links between points in two sets.
% One point in set A can only link to one point in set B or it can not be
% linked at all. The cost associated with the link from element i of A to
% element j of B is given by cc(i,j).
%
% SYNOPSIS [x, y] = lap(cc, NONLINK_MARKER, extendedTesting, augmentCC)
%
% INPUT:  cc: cost matrix, which has to be square. Set cost(i,j) to the
%             value of NONLINK_MARKER, if the link is not allowed. cc can
%             be input as a sparse matrix (in which case there won't be any
%             NONLINK_MARKERs).
%
%             For an unequal number of points (or, generally, if birth and
%             death is to be allowed) the cost matrix is formed as
%             a 2-by-2 catenation of four sub-matrices. If there are n
%             elements in set A and m elements in set B, sub-matrix (1,1)
%             is a n-by-m matrix with the cost for each link. The two
%             off-diagonal sub-matrices make non-links possible. They are
%             both diagonal square matrices (2,1: m-by-m; 1,2: n-by-n) with
%             the cost for not linking a point (e.g. determined by a
%             maximum search radius) in the diagonal and NONLINK_MARKER off
%             the diagonal. Sub-matrix (2,2) is a m-b-n matrix of any
%             not-NONLINK_MARKER value (suggested: zero, except for sparse
%             matrix input)
%
% NONLINK_MARKER : value to indicate that two points cannot be linked.
%             Default: -1. NaN is not allowed here
%
% extendedTesting: (optional, [0/{1}]) If 1, the cost matrix will be tested
%                  for violations (at the expense of speed):
%                     - There cannot be NaNs in cc
%                     - In every row and every column of cc there must be
%                       at least 1 element that is not a NONLINK_MARKER.
%
% augmentCC  (optional, [{0}/1]) If 1, the cost matrix will be augmented to
%            allow births and deaths.
%
% OUTPUT: x: The point A(i) links to B(x(i))
%         y: The point B(j) links to A(y(j))
%
%            Any x > m or y > n indicates that this point is not linked.
%
% lapjv is implemented through a MEX DLL.
%
% Author: Ge Yang, LCCB, TSRI, Dec. 10, 2004
% extended by jonas 6/05


% Basically, what this function does is to automatically extract the sparse matrix
% representation of cc.

% int *col = new int[size * NEIGHBOR_NUM_MAX + 1];
% int *first = new int[size + 2];
% int *x = new int[size + 1];
% int *y = new int[size + 1];
% double *u = new double[size + 1];
% double *v = new double[size + 1];
% double  *cc = new double[size * NEIGHBOR_NUM_MAX + 1];


%=====================
% TEST INPUT
%=====================
if (nargin == 1) || isempty(NONLINK_MARKER)
    NONLINK_MARKER = -1;
elseif isnan(NONLINK_MARKER)
    error('NONLINK_MARKER cannot be NaN!')
end
if nargin < 4 || isempty(augmentCC)
    augmentCC = 0;
end
% test size
scc = size(cc);
% check size only if no augmentation
if ~augmentCC
    if scc(1) ~= scc(2) || length(scc) > 2
        error('cost must be a 2D square matrix!')
    end
    
elseif length(scc) > 2
    error('cost must be a 2D matrix!')
else
    % if we're augmenting, sparse matrices are produced. This will be
    % problematic with cost 0
    if any(cc(:)==0)
        validCC = cc ~= NONLINK_MARKER;
        if any(any(cc(validCC))) < 0
            warning('there are negative costs. This could lead to errors')
        end
        cc(validCC) = cc(validCC) + 1;
    end
end

if nargin < 3 || isempty(extendedTesting)
    extendedTesting = 1;
end
%=======================

%=================================
% AUTMENT COST MATRIX IF SELECTED
%=================================
if augmentCC
    % expand the m-by-n cost matrix to a (m+n)-by-(n+m) matrix, adding
    % diagonals with a cost above the highest cost
    maxCost = max(cc(:)) + 1;

    % check if sparse
    if issparse(cc)
        % mmDiag = spdiags(maxCost * ones(scc(1),1), 0, scc(1), scc(1));
        % nnDiag = spdiags(maxCost * ones(scc(2),1), 0, scc(2), scc(2));
        % nmMat  = sparse(ones(scc(2), scc(1)));
        % cc = [cc, mmDiag; nnDiag, nmMat];
        cc = [cc, spdiags(maxCost * ones(scc(1),1), 0, scc(1), scc(1));...
            spdiags(maxCost * ones(scc(2),1), 0, scc(2), scc(2)),...
            sparse(ones(scc(2), scc(1)))];
    else
        % mmDiag = diag(maxCost * ones(scc(1),1));
        % nnDiag = diag(maxCost * ones(scc(2),1));
        % nmMat  = sparse(ones(scc(2), scc(1)));
        % cc = [cc, mmDiag; nnDiag, nmMat];
        
        % make cc sparse. Take NLM in cc into account!
        cc(cc==NONLINK_MARKER) = 0;
        cc = [cc, diag(maxCost * ones(scc(1),1)); ...
            diag(maxCost * ones(scc(2),1)), ...
            sparse(ones(scc(2), scc(1)))];
    end

    % remember that the size of the matrix has increased!
    scc = [sum(scc), sum(scc)];
end


%=============================
% CALCULATE SPARSE MATRIX
%=============================

% Calculate sparse representation of cost matrix with compactCC (all
% significant elements of cc), fst and kk. Read cc in rows, not cols!

% fst: for every first significant element of a row: index into compactCC
% kk : for every significant element: column

% do the work on the transposed cost matrix!
cc = cc';
% find the significant elements. If sparse input, find nonzero elements
if issparse(cc)
    [rowIdx, colIdx, val] = find(cc);
else
    [rowIdx, colIdx, val] = find(cc ~= NONLINK_MARKER);
end

% test that all cols and all rows are filled, and that there are no nans
if extendedTesting
    allCols = unique(colIdx);
    allRows = unique(rowIdx);
    if any(isnan(val)) || ~(length(allCols) == scc(1) && length(allRows) == scc(2))
        error('there must be at least one possible link per row and column, and there cannot be NaNs')
    end
end


% write value vector, pad a zero
compactCC = [0; val];
% write kk, pad a zero
kk = [0; rowIdx];

% colIdx is already sorted, so we can find out the number of entries per
% column via diff. Wherever there is a jump in the column, the deltaIdx
% will be 1, and its rowIdx will equal fst
deltaIdx = diff([0;colIdx]);
% add 0 and length+1
fst = [0;find(deltaIdx);length(val)+1];

%==================================


%==================================
% CALL MEX-FUNCTION
%==================================
[x, y, u, v] = mexLap(double(scc(1)), int32(length(compactCC)), double(compactCC), int32(kk), int32(fst));
%==================================

% remove first element from output vectors, as it is a meaningless 0.
x = x(2:end);
y = y(2:end);