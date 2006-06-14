
%for the initial simple linking from frame to frame
costMatrices(1).costMatFun = 'costMatSimple';
costMatrices(1).costMatParam = struct('searchRadius',2,'maxAmpRatio',2,'noLnkPrctl',-1);

%for linking between frames again using statistical data on the tracks
costMatrices(2).costMatFun = 'costMatLogL';
costMatrices(2).costMatParam = struct('cutoffProbD',0.05,'cutoffProbA',0.9999,'noLnkPrctl',-1);

%for gap closing
costMatrices(3).costMatFun = 'costMatCloseGaps';
costMatrices(3).costMatParam = struct('cutoffProbD1',0.15,'cutoffProbA1',0.9999,...
    'cutoffProbD2',0.1,'cutoffProbA2',0.99,'noLnkPrctl',-1);

%for merging and splitting
costMatrices(4).costMatFun = 'costVecLinkMS';
costMatrices(4).costMatParam = struct('cutoffProbD',0.15,'cutoffProbA',0.99);

%gap closing parameters
gapCloseParam.timeWindow = 10;
gapCloseParam.mergeSplit = 1;

%iteration parameters
iterParam.tolerance = 0.05;
iterParam.lenFrac = 0.3;

