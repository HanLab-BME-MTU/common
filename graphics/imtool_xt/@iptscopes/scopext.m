function scopext(ext)
%SCOPEXT  Register Image Processing Toolbox scope extensions.

%   Copyright 2007-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/04/25 07:11:27 $

r = ext.add('Tools',...
    'Image Tool',...
    'iptscopes.IMToolExport', ...
    getString(message('images:imtool:exportExtDesc')),...
    getString(message('images:commonUIString:imageTool')));
r.Depends = {'Visuals:Video'};

r = ext.add('Tools',...
    'Pixel Region',...
    'iptscopes.PixelRegion', ...
    getString(message('images:imtool:pixelRegionExtDesc')),...
    getString(message('images:commonUIString:pixelRegion')));
r.Depends = {'Visuals:Video'};

r = ext.add('Tools',...
    'Image Navigation Tools',...
    'iptscopes.IPTPanZoom', ...
    getString(message('images:imtool:panZoomExtDesc')),...
    getString(message('images:implayUIString:imageNavToolsLabel')));

r.Depends = {'Visuals:Video'};

% [EOF]