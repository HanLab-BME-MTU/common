function  = minWeightMatching(D, nonLinkMarker)

% D must be squared
% D must be symetric. Only the lower triangular part will be considered
% D can be either full or sparse
% D > 0

[n m] = size(D);

if n ~= m
  error('Input matrix must be squared.');
end

if any(D < 0)
  error('Weights must be >= 0');
end

[u, v, w] = find(tril(D));

M = minWeightMatchingMEX(n, [u v], w);

