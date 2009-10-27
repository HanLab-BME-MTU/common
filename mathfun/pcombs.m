function C = pcombs(v)
% This function generate all pairwise combinations of a vector v of length
% n. The number of possible pairwise combinations is n (n - 1) / 2. The
% elements of v must be all different.
%
% e.g.:
% v = 1 2 4 3
% C = 
%   1  2
%   1  4
%   2  4
%   1  3
%   2  3
%   4  3

assert(min(v == unique(v)) == 1);

[I J] = meshgrid(v, v');
I = triu(I, 1);
J = triu(J, 1);
C = horzcat(nonzeros(J), nonzeros(I));

end