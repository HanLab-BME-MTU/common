function [x, y] = lap(cc, NONLINK_MARKER, extendedTesting)
%LAP solves the linear assignment problem for a given cost matrix
%
% A linear assignment tries to establish links between points in two sets.
% One point in set A can only link to one point in set B or it can not be
% linked at all. The cost associated with the link from element i of A to
% element j of B is given by cc(i,j).
%
% SYNOPSIS [x, y] = lap(cc, NONLINK_MARKER, extendedTesting)
%
% INPUT:  cc: cost matrix, which has to be square. Set cost(i,j) to the
%             value of NONLINK_MARKER, if the link is not allowed. %
%
%             For an unequal number of points the cost matrix is formed as
%             a 2-by-2 catenation of four sub-matrices. If there are n
%             elements in set A and m elements in set B, sub-matrix (1,1)
%             is a n-by-m matrix with the cost for each link. The two
%             off-diagonal sub-matrices make non-links possible. They are
%             both diagonal square matrices (2,1: m-by-m; 1,2: n-by-n) with
%             the cost for not linking a point (e.g. determined by a
%             maximum search radius) in the diagonal and NONLINK_MARKER off
%             the diagonal. Sub-matrix (2,2) is a m-b-n matrix of zeros.
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
% OUTPUT: x: The point A(i) links to B(x(i))
%         y: The point B(j) links to A(y(j))
%
%            Any x > m or y > n indicates that this point is not linked.
%
% lapjv is implemented through a MEX DLL.
%
% Author: Ge Yang, LCCB, TSRI, Dec. 10, 2004


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
if (nargin == 1)
    NONLINK_MARKER = -1;
elseif isnan(NONLINK_MARKER)
    error('NONLINK_MARKER cannot be NaN!')
end


scc = size(cc);
if scc(1) ~= scc(2) || length(scc) > 2
    error('cost must be a 2D square matrixt!')
end

if nargin < 3 || isempty(extendedTesting)
    extendedTesting = 1;
end
%=======================


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
    [rowIdx, colIdx, val] = find(cc ~= NONLINK_MARKER),
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