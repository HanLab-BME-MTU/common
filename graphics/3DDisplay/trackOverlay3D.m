function [overlays]=trackOverlay3D(MD,tracksFinal,processedFrames)
tracksColors=uint8(hsv*255);

overlays=cell(1,numel(processedFrames));
[xSize,ySize,zSize]=size(MD.getChannel(1).loadStack(1));

parfor frameIdx=1:numel(overlays)
    disp(['Generate overlay for frame ' int2str(frameIdx)])
    overlay=zeros(xSize,ySize,3,zSize);
    for trackIdx=1:length(tracksFinal)
        trackFrames=getTrackFrames(tracksFinal,trackIdx);
        if(ismember(frameIdx,trackFrames))
            RGB=tracksColors(mod(trackIdx,length(tracksColors))+1,:);
            yCoord=tracksFinal(trackIdx).tracksCoordAmpCG(:,1:8:(frameIdx-trackFrames+1)*8); 
            xCoord=tracksFinal(trackIdx).tracksCoordAmpCG(:,2:8:(frameIdx-trackFrames+1)*8); 
            zCoord=tracksFinal(trackIdx).tracksCoordAmpCG(:,3:8:(frameIdx-trackFrames+1)*8); 
            sampling=100;
            xSeg=max(1,round(linspaceNDim(xCoord(1:(end-1)),xCoord(2:(end)),sampling)));
            ySeg=max(1,round(linspaceNDim(yCoord(1:(end-1)),yCoord(2:(end)),sampling)));
            zSeg=max(1,round(linspaceNDim(zCoord(1:(end-1)),zCoord(2:(end)),sampling)));
            indx=sub2ind(size(overlay),xSeg,ySeg,ones(size(xSeg)),zSeg);
            overlay(indx)=RGB(1);
            indx=sub2ind(size(overlay),xSeg,ySeg,2*ones(size(xSeg)),zSeg);
            overlay(indx)=RGB(2);
            indx=sub2ind(size(overlay),xSeg,ySeg,3*ones(size(xSeg)),zSeg);
            overlay(indx)=RGB(3);
        end
    end
    overlays{frameIdx}=overlay;
end



function [trackFrames] = getTrackFrames (tracksFinal,trackIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;

trackFrames=(startTime:endTime)';
