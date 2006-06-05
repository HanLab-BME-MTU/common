
%for the initial simple linking from frame to frame
costMatrices(1).costMatFun = 'costMatSimple';
costMatrices(1).costMatParam = struct('searchRadius',3,'maxAmpRatio',2,'noLnkPrctl',-1);

%for linking between frames again using statistical data on the tracks
costMatrices(2).costMatFun = 'costMatLogL';
costMatrices(2).costMatParam = struct('cutoffCProb',0.9999,'noLnkPrctl',-1);

%for gap closing
costMatrices(3).costMatFun = 'costMatCloseGaps';
costMatrices(3).costMatParam = struct('cutCProb1',0.9999,'cutCProb2',0.999,'noLnkPrctl',-1);

%for merging and splitting
costMatrices(4).costMatFun = 'costVecLinkMS';
costMatrices(4).costMatParam = struct('cutoffCProb',0.99);

