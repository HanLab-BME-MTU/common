function [ max_v, max_i ] = ppmax( pp )
%ppmax Find the absolute maximum of a piecewise polynomial
%
% INPUT
% pp - piecewise polynomial
%
% OUTPUT
% max_v - value of maximum on the interval where pp is defined
% max_i - abscissa where absolute maximum occurs
%
% See also mkpp, ppdiff, pproots, ppextrema
%
% Mark Kittisopikul, Februrary 2017
% Jaqaman Lab
% UT Southwestern

maxima = ppextrema(pp);


maxima = [pp.breaks(1); maxima ; pp.breaks(end)];
[max_v, max_i] = max(ppval(pp, maxima));
max_i = maxima(max_i);


end

