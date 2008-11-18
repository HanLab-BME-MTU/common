function trueOrFalse = isCalledFromBase
%ISCALLEDFROMBASE checks whether a function has been called from the command line directly
%
% SYNOPSIS: trueOrFalse = isCalledFromBase
%
% INPUT none
%
% OUTPUT trueOrFalse: true if the function's caller is the base workspace
%
% REMARKS
%
% created with MATLAB ver.: 7.7.0.471 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 17-Nov-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a variable in the calling workspace
assignin('caller','callerCheck1234567890',1);
% try to find that variable in the base workspace
trueOrFalse = evalin('base','exist(''callerCheck'',''var'');');
% remove the variable from the calling workspace
evalin('caller','clear callerCheck1234567890')