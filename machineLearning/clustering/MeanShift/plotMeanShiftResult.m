function fh = plotMeanShiftResult(data,clusterInfo,ptMap,varargin)
%PLOTMEANSHIFTRESULT displays results of mean shift clustering
%   
%   Input:
%              data -> data to be clustered
%       clusterInfo -> cluster information from mean shift clustering
%             ptMap -> point-to-cluster mapping from mean shift clustering
%
%   Output:
%       fh -> figure handle
%
%   US, 2012/11/09
%

ip=inputParser;

ip.addRequired('data',@isnumeric);
ip.addRequired('clusterInfo',@isstruct);
ip.addRequired('ptMap',@isnumeric);

ip.parse(data,clusterInfo,ptMap,varargin{:});

data=ip.Results.data;
clusterInfo=ip.Results.clusterInfo;
ptMap=ip.Results.ptMap;

fh=figure;
hold on;

clusterColors=mat2gray(squeeze(label2rgb(1:numel(clusterInfo))),[0,255]);
    
for k = 1:numel(clusterInfo)
    ptCurClusterCenter=clusterInfo(k).ptClusterCenter;
    plot(data(ptMap==k,1),data(ptMap==k,2), ...
          '.', 'Color', clusterColors(k,:))
    plot(ptCurClusterCenter(1),ptCurClusterCenter(2), ...
        'o','MarkerEdgeColor','k','MarkerFaceColor',clusterColors(k,:), 'MarkerSize',10)
end

hold off;

end

