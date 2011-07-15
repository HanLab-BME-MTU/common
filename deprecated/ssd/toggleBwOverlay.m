function toggleBwOverlay()
%TOGGLEBWOVERLAY toggles the B & W overlay status
%
% SYNOPSIS toggleBwOverlay
%
% INPUT none;
%
% OUTPUT none;

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);


global bwOverlays__;

bwOverlays__ = mod(bwOverlays__ + 1,2);