function [overlays]=featureOverlay3D(MD,tracksFinal,processedFrames)
% IN DEV. DONT COMMIT

overlays=cell(1,numel(processedFrames));
[xSize,ySize,zSize]=size(MD.getChannel(1).loadStack(1))

for frameIdx=1:numel(overlays)
    disp(['Printing frame ' int2str(frameIdx)])
    overlay=zeros(xSize,ySize,3,zSize);
    for trackIdx=1:length(tracksFinal)
        trackFrames=getTrackFrames(tracksFinal,trackIdx);
        if(ismember(frameIdx,trackFrames))
            disp(['Printing track ' int2str(trackIdx)])
            yCoord=tracksFinal(trackIdx).tracksCoordAmpCG(:,1:8:(frameIdx-trackFrames+1)*8); 
            xCoord=tracksFinal(trackIdx).tracksCoordAmpCG(:,2:8:(frameIdx-trackFrames+1)*8); 
            zCoord=tracksFinal(trackIdx).tracksCoordAmpCG(:,3:8:(frameIdx-trackFrames+1)*8); 
            indx=sub2ind(size(overlay), ...
            [max(1,floor(xCoord)) max(1,floor(xCoord)-1) ], ...print
            [max(1,floor(yCoord)) max(1,floor(yCoord)-1) ], ...A
            [ones(size(xCoord)) ones(size(xCoord))],...
            [max(1,floor(zCoord)) max(1,floor(zCoord)-1) ]...
            );
            overlay(indx)=255;
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
