function printOverlay3D(MD,overlays,folder,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('saturation',[0 1],@isnumeric);
ip.parse(varargin{:});

parfor i=1:numel(overlays)
    vol=MD.getChannel(1).loadStack(i);
    volRGB=repmat(uint8(255*mat2gray(vol,quantile(double(vol(:)),ip.Results.saturation))),[1 1 1 3]);
    volRGB=permute(volRGB,[1 2 4 3]);
    volRGB(overlays{i}>0)=overlays{i}(overlays{i}>0);
    fname=[folder filesep  '3DOverlays' num2str(i,'%04.f') '.tif'];
    stackWrite(volRGB,fname)
end
    