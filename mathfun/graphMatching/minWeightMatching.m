function [U V] = minWeightMatching(D, nonLinkMarker)

% D must be squared
% D must be symetric. Only the lower triangular part will be considered
% D can be either full or sparse
% D > 0

[n m] = size(D);

if nargin < 2 || isemtpy(nonLinkMarker)
  nonLinkMarker = max(D(:)) + 1;
end

if n ~= m
  error('Input matrix must be squared.');
end

if any(D < 0)
  error('Weights must be >= 0');
end

% FIXME

[u, v, w] = find(tril(D));

M = perfectMatchingMEX(n, [u v], w);

% TODO: remove edge where w == nonLinkMarker