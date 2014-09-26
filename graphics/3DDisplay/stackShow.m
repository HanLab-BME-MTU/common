function stackShow(vol,varargin)
% Summary: Plot 3D stack <vol> as a set of voxels with optional overlay.
%
% Notes:
%   - wrapper for vol3D ((C) Conti, Woodford)
%   - Voxel values are used for alpha.  
%   - Probably faster than viewStack (rough test)
%   - Allow the use of an overlay mask in Red 
%   
% Input:
%   - vol: 3D matrix 
%
% Options:
%   - 'overlay': overlay matrix (such as mask) same size of vol (NOT
%   ASSESSED). Non-zero value results in a Red mask. 
%   - 'threshold': optional hard treshold on voxel plotting
%
% Output:
%   - None
%
% Philippe Roudot 09/2014

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol',@isnumeric);
ip.addParamValue('overlay',[],@isnumeric);
ip.addParamValue('threshold',0,@isnumeric);
ip.parse(vol, varargin{:});

vol(vol<ip.Results.threshold)=0;

if ~isempty(ip.Results.overlay)
    label=ip.Results.overlay;
    vol=vol/max(vol(:));
    volView=repmat(vol,[1 1 1 3]);
    volView=double(volView);
    volView(:,:,:,1)=label;
    vol3d('cdata',volView,'Alpha',vol);     
else
    vol=vol/max(vol(:));
    volView=repmat(vol,[1 1 1 3]);
    volView=double(volView);
    vol3d('cdata',volView,'Alpha',vol);   
    colormap('gray');
end    

% Respect proportional aspect ratios
daspect([1 1 1]);                      