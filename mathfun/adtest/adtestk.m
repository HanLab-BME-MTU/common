%[hval T cval A2] = adtestk(samples)
%
% Inputs:
%         samples : cell array of sample vectors (arrays)
%
% Outputs:
%            hval : Result of the hypothesis test
%                   0: Do not reject H0
%                   1: Reject H0
%               T : value of the test statistic
%            cval : critical value
%              A2 : value of the k-sample distance

% Francois Aguet, 03/08/2012

function [hval, T, cval, A2] = adtestk(samples) %#ok<STOUT,INUSD>