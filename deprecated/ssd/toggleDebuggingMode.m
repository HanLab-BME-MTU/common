function toggleDebuggingMode()
%TOGGLEDEBUGGINGMODE toggles the debuggingMode status for G.D. functions
%
% SYNOPSIS toggleDebuggingMode
%
% INPUT none;
%
% OUTPUT none;

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

global debuggingMode__;

debuggingMode__ = mod(debuggingMode__ + 1,2);