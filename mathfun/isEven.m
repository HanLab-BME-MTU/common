function out = isEven(in)
%ISEVEN checks whether a number is even
%
% SYNOPSIS out = isEven(in)
%
% INPUT    in :  input (array) of numbers to be tested. IsEven only works properly for integers. 
% OUTPUT   out:  array of size(in) with 
%                   1 for even integers and zero
%                   0 for odd integers
%                 NaN for non-integers
%
% c: jonas 5/05
% Last modified 10/16/2009 - Francois Aguet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NaNs cannont be converted to logical type (true or false only)
out = mod(in+1, 2);
out((out ~= 0) & (out ~= 1)) = NaN;