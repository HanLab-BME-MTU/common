function [areaCell,areaConvert] = getCellArea(movieData)
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

if sum(maskProc.checkChannelOutput)==1
    iMaskChan = find(maskProc.checkChannelOutput);
    maskCell = maskProc.loadChannelOutput(iMaskChan,1);
else
    I=double(movieData.channels_(1).loadImage(1));
    %Combine the the multiple masks to one
    maskEach = arrayfun(@(x) maskProc.loadChannelOutput(x,1),find(maskProc.checkChannelOutput),'UniformOutput',false);
%     maskAll=reshape(cell2mat(maskEach),size(I,1),size(I,2),[]);
    areaEach = cellfun(@(x)sum(x(:)),maskEach);
    [~,iMin] = min(areaEach);
    maskCell = maskEach{iMin};
%     maskCell = all(maskAll,3); %changed from any
end

pixSize_mu=movieData.pixelSize_*1e-3; % in um/pixel
areaConvert=pixSize_mu^2; % in um2/pixel

areaCell = sum(maskCell(:))*areaConvert;  % in um2
