function [x, y] = lap(cc, NONLINK_MARKER)
%LAP solves the linear assignment problem for a given cost matrix
%
% A linear assignment tries to establish links between points in two sets.
% One point in set A can only link to one point in set B or it can not be
% linked at all. The cost associated with the link from element i of A to
% element j of B is given by cc(i,j). 
%
% SYNOPSIS [x, y] = lap(cc, NONLINK_MARKER)
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

%=======================


%=============================
% CALCULATE SPARSE MATRIX
%=============================

% Calculate sparse representation of cost matrix with compactCC (all
% significant elements of cc), fst and kk. Read cc in rows, not cols!

% fst: for every first significant element of a row: index into compactCC
% kk : for every significant element: column

% logical matrix with 1 for significant elements
cc = cc';
logicArrayT = (cc ~= NONLINK_MARKER); % transpose here?
logicArray = logicArrayT';

rowSumLA = sum(logicArray,2);
colSumLA = sum(logicArray,1);
if any(isnan(rowSumLA) | rowSumLA == 0) || any(isnan(colSumLA) | colSumLA == 0)
    error('there must be at least one possible link per row and column, and there cannot be NaNs')
end

% assign all significant elements
compactCC_2 = [0;cc(logicArrayT)];
% find all column indices, pad a zero. 
[kk_2,dummy] = find(logicArrayT);
kk_2 = [0;kk_2];
% find all the first entries per row
fst_2 = cumsum(rowSumLA); 
fst_2 = [-1;0;fst_2] + 1; 

%==================================


%==================================
% CALL MEX-FUNCTION
%==================================
[x, y, u, v] = mexLap(double(scc(1)), int32(length(compactCC_2)), double(compactCC_2), int32(kk_2), int32(fst_2));
%==================================

% remove first element from output vectors, as it is a meaningless 0.
x = x(2:end);
y = y(2:end);