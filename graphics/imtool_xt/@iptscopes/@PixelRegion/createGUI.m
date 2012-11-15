function plugInGUI = CreateGUI(this)   %#ok hSrc not used
%CreateGUI Build and cache UI plug-in for PixelRegion plug-in.
%   This adds the button and menu to the scope.
%   No install/render needs to be done here.

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/08/09 17:55:18 $

hSrc = this.Application.DataSource;

if isempty(hSrc) || ~isDataLoaded(hSrc)
    enab = 'off';
else
    enab = 'on';
end

mPixel = uimgr.uimenu('PixelRegion', 1, getString(message('images:imtoolUIString:pixelRegionMenubarLabel')));
mPixel.Enable = enab;
mPixel.setWidgetPropertyDefault(...
    'busyaction', 'cancel', ...
    'callback', @(hco, ev) launch(this));

bPixel = uimgr.uipushtool('PixelRegion', 1);
bPixel.IconAppData = 'pixel_region';
bPixel.Enable = enab;
bPixel.setWidgetPropertyDefault(...
    'busyaction', 'cancel', ...
    'tooltip', getString(message('images:imtoolUIString:pixelRegionTooltipString')), ...
    'click',   @(hco, ev) launch(this));

% Create plug-in installer
plan = {mPixel, 'Base/Menus/Tools/VideoTools'; ...
        bPixel, 'Base/Toolbars/Main/Tools/Standard'};
plugInGUI = uimgr.Installer(plan);

% [EOF]
