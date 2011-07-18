function y = clipNegative(x)
%CLIPNEGATIVE clips all negative values to zero
%
% SYNOPSIS y = clipNegative(x)
%
%    where y = 0 if x < 0; y = x otherwise

fprintf(2,['Warning: ''' mfilename ''' is deprecated. Use normpdf instead.\n']);

aux = (x >= 0);
y = aux .* x;
