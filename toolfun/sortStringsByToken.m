%[s, rm] = sortStringsByToken(s, token, mode) sorts input array by a number sequence preceding or following a specified token
%
% Inputs:
%      s : cell array of strings    
%  token : token to match
%   mode : 'pre' (sort by number sequence preceding token) or
%          'post' (by sequence following token)
%
% Outputs:
%      s : sorted array
%     rm : index of input that did not match token and was excluded from output

% Francois Aguet, 02/08/2013

function [s, rm] = sortStringsByToken(s, token, mode)

switch mode
    case 'post'
        idx = cellfun(@(x) str2double(regexpi(x,['(?<=' token ')\d+'], 'match')), s, 'UniformOutput', false);
    case 'pre'
        idx = cellfun(@(x) str2double(regexpi(x,['\d+(?=' token ')'], 'match')), s, 'UniformOutput', false);
end
rm = find(cellfun(@isempty, idx));
s(rm) = [];
[~,idx] = sort([idx{:}]);
if ~isempty(idx)
    s = s(idx);
end