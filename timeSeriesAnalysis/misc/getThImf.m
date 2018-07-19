function [out] = getThImf(imf,k)
%Stupid function to get ith imf from cell array
% That is because matlab does not allow if in an anomimous function
%Marco Vilela, 2012
out = [];
if size(imf,1) > k
    out = imf(k,:);
end