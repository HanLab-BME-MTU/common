function [stack, stackinfo] = stackRead(stackpath)
%
% [stack, stackinfo] = stackRead(stackpath)
%
% Wrapper function for Fran�ois N�d�lec's tiffread.m
% This works on STK files as well as multipage TIFFs.
% Outputs:
%
%   stack     : 3D data array
%   stackinfo : Any other information included in original file
%
% Fran�ois Aguet, 01/2010

stackinfo = tiffread(stackpath);
stack = cat(3,stackinfo(:).data);
stackinfo = rmfield(stackinfo, 'data');