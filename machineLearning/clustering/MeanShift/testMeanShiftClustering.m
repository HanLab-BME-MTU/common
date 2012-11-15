clc
clear
close all

% generate three synthetic gaussian clusters
numPointsPerCluster = 2^16;
ptTrueClusterCenters = [ 1 1 ; -1, -1; 1, -1 ];
clusterStdDev = 0.6;

ptRawData = [];
for i = 1:size(ptTrueClusterCenters,1)
    
    ptCurCluster = repmat(ptTrueClusterCenters(i,:), numPointsPerCluster, 1) + ...
                   clusterStdDev * randn(numPointsPerCluster,2);               
    ptRawData = [ ptRawData; ptCurCluster ];
    
end

% Run mean-shift
bandwidth = 0.5;
tic
profile on;
[clusterInfo,pointToClusterMap] = MeanShiftClustering(ptRawData, bandwidth, ... 
                                                      'flagDebug', false, ...
                                                      'kernel', 'gaussian', ...
                                                      'flagUseKDTree', true);
profile off;
profile viewer;
toc

% plot result
figure;
hold on;

    clusterColors = mat2gray( squeeze( label2rgb( 1:numel(clusterInfo) ) ), [0,255] );
    for k = 1:numel(clusterInfo)
        ptCurClusterCenter = clusterInfo(k).ptClusterCenter;
        plot( ptRawData(pointToClusterMap==k, 1), ...
              ptRawData(pointToClusterMap==k, 2), ...
              '.', 'Color', clusterColors(k,:))
        plot(ptCurClusterCenter(1),ptCurClusterCenter(2),'o','MarkerEdgeColor','k','MarkerFaceColor',clusterColors(k,:), 'MarkerSize',10)
    end
    title( sprintf('Mean-shift clustering result - %d clusters were found', numel(clusterInfo)));

hold off;
