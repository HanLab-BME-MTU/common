
function runBatchDetectTrack(batchJob)

%get number of jobs
numJobs = length(batchJob);

for iJob = 1 : numJobs
    
    disp(['--- Movie #' num2str(iJob) ' ---'])
    
    %detect
    try    
    [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
        detectSubResFeatures2D_StandAlone(batchJob(iJob).movieParam,...
        batchJob(iJob).detectionParam,batchJob(iJob).saveResultsDetect);
    clear background exceptions localMaxima psfSigma 
    catch
        disp(['Job ' num2str(iJob) ': Detection failure']);
    end
    
    %track
    try
    [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(1:200),...
        batchJob(iJob).costMatrices,batchJob(iJob).gapCloseParam,...
        batchJob(iJob).kalmanFunctions,batchJob(iJob).probDim,...
        batchJob(iJob).saveResultsTrack,batchJob(iJob).verbose);
    clear tracksFinal kalmanInfoLink errFlag movieInfo
    catch
        disp(['Job ' num2str(iJob) ':Tracking failure']);
        clear movieInfo
    end
    
end

