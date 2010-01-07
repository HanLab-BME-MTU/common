function metaData = makeMetaData(sampArray,nReinit)

%This function takes all the protrusion and activity samples from the movies stored in
%sampArray and concatenates them into one big dataset

if nargin < 2
    nReinit = [];
end

nMovies = length(sampArray);

%Find maximum number of bands, frames
for j = 1:nMovies
    [movNband(j),movNwin(j),movNframes(j)] = size(sampArray(j).activity.avgActivity);
end

nBands = max(movNband);

nFramesMax = max(movNframes);

metaData.avgActivity = nan(nBands,0,nFramesMax);
metaData.nPixels = nan(nBands,0,nFramesMax);
metaData.stdActivity = nan(nBands,0,nFramesMax);
metaData.meanPos = nan(nBands,0,nFramesMax,2);
metaData.nPoints = nan(0,nFramesMax);
metaData.averageNormalComponent = nan(0,nFramesMax);
metaData.averageMagnitude = nan(0,nFramesMax);
metaData.stdNormalComponent = nan(0,nFramesMax);
metaData.stdMagnitude = nan(0,nFramesMax);

%Concatenate all the windows
for iMov = 1:nMovies
    
    %Concatenate activity samples
    tmpDat = cat(3,sampArray(iMov).activity.avgActivity,nan(movNband(iMov),movNwin(iMov),nFramesMax-movNframes(iMov)));
    tmpDat = cat(1,tmpDat,nan(nBands-movNband(iMov),movNwin(iMov),nFramesMax));    
    metaData.avgActivity = cat(2,metaData.avgActivity,tmpDat);    
    
    tmpDat = cat(3,sampArray(iMov).activity.nPixels,nan(movNband(iMov),movNwin(iMov),nFramesMax-movNframes(iMov)));
    tmpDat = cat(1,tmpDat,nan(nBands-movNband(iMov),movNwin(iMov),nFramesMax));
    metaData.nPixels = cat(2,metaData.nPixels,tmpDat);    
    
    tmpDat = cat(3,sampArray(iMov).activity.stdActivity,nan(movNband(iMov),movNwin(iMov),nFramesMax-movNframes(iMov)));
    tmpDat = cat(1,tmpDat,nan(nBands-movNband(iMov),movNwin(iMov),nFramesMax));
    metaData.stdActivity = cat(2,metaData.stdActivity,tmpDat);    
    
    tmpDat = cat(3,sampArray(iMov).activity.meanPos,nan(movNband(iMov),movNwin(iMov),nFramesMax-movNframes(iMov),2));
    tmpDat = cat(1,tmpDat,nan(nBands-movNband(iMov),movNwin(iMov),nFramesMax,2));
    metaData.meanPos = cat(2,metaData.meanPos,tmpDat);  
    
    %Concatenate protrusion samples
    tmpDat = cat(2,sampArray(iMov).protrusion.nPoints,nan(movNwin(iMov),nFramesMax-movNframes(iMov)+1));
    metaData.nPoints = cat(1,metaData.nPoints,tmpDat);
    
    tmpDat = cat(2,sampArray(iMov).protrusion.averageNormalComponent,nan(movNwin(iMov),nFramesMax-movNframes(iMov)+1));
    metaData.averageNormalComponent = cat(1,metaData.averageNormalComponent,tmpDat);
    
    tmpDat = cat(2,sampArray(iMov).protrusion.averageMagnitude,nan(movNwin(iMov),nFramesMax-movNframes(iMov)+1));
    metaData.averageMagnitude = cat(1,metaData.averageMagnitude,tmpDat);
    
    tmpDat = cat(2,sampArray(iMov).protrusion.stdNormalComponent,nan(movNwin(iMov),nFramesMax-movNframes(iMov)+1));
    metaData.stdNormalComponent = cat(1,metaData.stdNormalComponent,tmpDat);
    
    tmpDat = cat(2,sampArray(iMov).protrusion.stdMagnitude,nan(movNwin(iMov),nFramesMax-movNframes(iMov)+1));
    metaData.stdMagnitude = cat(1,metaData.stdMagnitude,tmpDat);
    
    
end

%Seperate windows after re-initialization
if ~isempty(nReinit) && (nReinit < min(movNframes))
    
    newMeta.avgActivity = nan(nBands,0,nReinit);
    newMeta.nPixels = nan(nBands,0,nReinit);
    newMeta.stdActivity = nan(nBands,0,nReinit);
    newMeta.meanPos = nan(nBands,0,nReinit,2);
    newMeta.nPoints = nan(0,nReinit);    
    newMeta.averageNormalComponent = nan(0,nReinit);
    newMeta.averageMagnitude = nan(0,nReinit);
    newMeta.stdNormalComponent = nan(0,nReinit);
    newMeta.stdMagnitude = nan(0,nReinit);
    
    
    iReinit = 0:nReinit:nFramesMax;

    
     if iReinit(end) ~= nFramesMax
         %iReinit = [iReinit nFramesMax];
         disp(['Warning: ' num2str(nFramesMax-iReinit(end)) ' frames dropped to make total length a multiple of nReinit!']);
     end
    
    for j = 1:(length(iReinit)-1)
        
        %Activity
        newMeta.avgActivity = cat(2,newMeta.avgActivity,metaData.avgActivity(:,:,iReinit(j)+1:iReinit(j+1)));
        newMeta.nPixels = cat(2,newMeta.nPixels,metaData.nPixels(:,:,iReinit(j)+1:iReinit(j+1)));   
        newMeta.stdActivity = cat(2,newMeta.stdActivity,metaData.stdActivity(:,:,iReinit(j)+1:iReinit(j+1)));
        newMeta.meanPos = cat(2,newMeta.meanPos,metaData.meanPos(:,:,iReinit(j)+1:iReinit(j+1),:));
        
        %Protrusion
        newMeta.nPoints = cat(1,newMeta.nPoints,metaData.nPoints(:,iReinit(j)+1:iReinit(j+1)));
        newMeta.averageNormalComponent = cat(1,newMeta.averageNormalComponent,metaData.averageNormalComponent(:,iReinit(j)+1:iReinit(j+1)));       
        newMeta.averageMagnitude = cat(1,newMeta.averageMagnitude,metaData.averageMagnitude(:,iReinit(j)+1:iReinit(j+1)));
        newMeta.stdNormalComponent = cat(1,newMeta.stdNormalComponent,metaData.stdNormalComponent(:,iReinit(j)+1:iReinit(j+1)));
        newMeta.stdMagnitude = cat(1,newMeta.stdMagnitude,metaData.stdMagnitude(:,iReinit(j)+1:iReinit(j+1)));
        
    end

    metaData.avgActivity = newMeta.avgActivity;
    metaData.nPixels = newMeta.nPixels;
    metaData.stdActivity = newMeta.stdActivity;
    metaData.meanPos = newMeta.meanPos;
    metaData.nPoints = newMeta.nPoints;
    metaData.averageNormalComponent = newMeta.averageNormalComponent;
    metaData.averageMagnitude = newMeta.averageMagnitude;
    metaData.stdNormalComponent = newMeta.stdNormalComponent;
    metaData.stdMagnitude = newMeta.stdMagnitude;
    
end

