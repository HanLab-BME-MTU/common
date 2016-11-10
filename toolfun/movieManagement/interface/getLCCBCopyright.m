function copyright = getLCCBCopyright()
%
% This is a user-defined function used in UTSW software. 
% It is called when any GUI is generated. It configures the copyright
% information.
%
% Input: 
%
%
% Output:
%
%   copyright - String: copyright and version information
%
% Chuangang Ren, 11/2010
% Sebastien Besson, Feb 2013
% Andrew Jamieson, Nov 2016 - UTSW

% Set year and version information
str_year = '2016';
copyright = sprintf('Copyright %s UTSW', str_year);
