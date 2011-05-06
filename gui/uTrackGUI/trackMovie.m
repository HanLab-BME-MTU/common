function trackMovie(movieData,paramsIn)
% Track features in a movie which has been processed by a detection method
%
% Sebastien Besson, 5/2011

%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('TrackingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(TrackingProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

obj = movieData.processes_{iProc};

iDetection =movieData.getProcessIndex('DetectionProcess',1,0);
detectionProcess = movieData.processes_{iDetection};
detectionProcess.checkChannelOutput(trackingProcess.channelIndex_);

obj.setInFilePaths(detectionProcess.outFilePaths_);

for i = obj.channelIndex_
    
    load(obj.inFilePaths_{i},'movieInfo');
    obj.funParams_.saveResults.filename = ['Channel_' num2str(i) '_' obj.filename_];
    
    %Check/create directory
    if ~exist(obj.funParams_.saveResults.dir,'dir')
        mkdir(obj.funParams_.saveResults.dir)
    end
    
    if ~obj.overwrite_
        % file name enumeration
        obj.funParams_.saveResults.filename = enumFileName(obj.funParams_.saveResults.dir, obj.funParams_.saveResults.filename);
    end
    
    % Call function - return tracksFinal for reuse in the export
    % feature
    tracksFinal = trackCloseGapsKalmanSparse(movieInfo, obj.funParams_.costMatrices, obj.funParams_.gapCloseParam,...
        obj.funParams_.kalmanFunctions, obj.funParams_.probDim, obj.funParams_.saveResults, obj.funParams_.verbose);
    
    obj.setOutFilePath(i,[obj.funParams_.saveResults.dir filesep obj.funParams_.saveResults.filename]);
    
    % Optional export
    if obj.funParams_.saveResults.export
        if ~obj.funParams_.gapCloseParam.mergeSplit
            [M.trackedFeatureInfo M.trackedFeatureIndx]=...
                convStruct2MatNoMS(tracksFinal);
        else
            [M.trackedFeatureInfo M.trackedFeatureIndx,M.trackStartRow,M.numSegments]=...
                convStruct2MatIgnoreMS(tracksFinal);
        end
        
        matResultsSaveFile=[obj.funParams_.saveResults.dir filesep obj.funParams_.saveResults.filename(1:end-4) '_mat.mat'];
        save(matResultsSaveFile,'-struct','M');
        clear M;
    end
end

obj.setDateTime;
movieData.save; %Save the new movieData to disk
end