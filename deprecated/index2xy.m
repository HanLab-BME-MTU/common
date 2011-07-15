function [x y]=index2xy(ind,L,W)

% DESCRIPTION: convert from pixel index to xy coordinates
%
% SYNOPSIS: [x y]=index2xy(ind,L,W)
%
% INPUT: ind: vector of indices to be converted
%        L: image length
%        W: image width
%
% OUTPUT: x: vector of x (column) coordinates
%         y: vector of y (row) coordinates
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
% USERNAME: kathomps
% DATE: 31-Oct-2006
%
%

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

x=ceil(ind./L);
y=ind-(x-1)*L;
