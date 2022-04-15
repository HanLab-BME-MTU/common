function areaCell = getCellArea(movieData)
% function area = getCellArea(MD) returns area from mask refinement process
% of MD
% input:
%       movieData       movieData
% output:
%       areaCell        area in um2
% Sangyoon Han, Apr 2022

iMaskProcess=movieData.getProcessIndex('MaskRefinementProcess');
if ~isempty(iMaskProcess)
    maskProc = movieData.getProcess(iMaskProcess);
else
    error('Please run Mask Refinement Process and run this function again')
end
maskCell = maskProc.loadChannelOutput(iChan,ii);

pixSize_mu=movieData.pixelSize_*1e-3; % in um/pixel
areaConvert=pixSize_mu^2; % in um2/pixel

areaCell = sum(maskCell(:))*areaConvert;  % in um2
