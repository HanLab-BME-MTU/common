function [clusterInfo, pointToClusterMap, pointTraj] = VariableBandwidthMeanShiftClustering( ptData, bandwidth, varargin )
% Implementation of the variable bandwidth mean-shift clustering algorithm using KD-Trees
% 
%     Required Input Arguments:
% 
%                     ptData: specifies the input pointset to be clustered
%                             must be a numPoints x numDimension matrix 
% 
%                  bandwidth: size of the bandwidth to be used in the mean-shift kernels
%                             It must be either a a numPoints x 1 vector or
%                             a numDimension x numDimension x numPoints
%                             matrix
%                             
%     Optional Input Arguments:
%         
%                     kernel: specifies the type of kernel used in the non-parametric 
%                             (parzen-window) density function                 
%                              Default: 'gaussian'
%                              Options:                  
%                                 'gaussian' -- uses a gaussian kernel
%                                 'flat' - uses a flat box averaging kernel
%                     
%                     method: specifies the type of mean-shift algorithm to use
%                             Default: 'optimized'
%                             Options:
%                                  'standard': mode-seeking is done for each and every point separately
%                                 'optimized': *NOT YET SUPPORTED!* all nearby points encountered in the mode-seeking path of a point
%                                              are recorded as visited and mode-seeking is not done for them 
%                                              separately. This gives rise to a significant amount of speedup.
%                                  
%                                  
%              maxIterations: specifies the maximum number of iterations/mean-shifts 
%                             that are allowed searching for the mode of a point                        
%                             Default: 500
%         
%         minClusterDistance: specifies the minimum "distance" between two
%                             cluster centers before they are considered a
%                             single cluster. This distance is actually a
%                             percentile using an elliptical gaussian
%                             approximation based on mahalanobis distance.
%
%              kernelSupport: specifies the maximum distance of points to
%                             include when evaluating the mean-shift kernel.
%                             Default: =2*bandwidth for 'gaussian' and user-specified kernel
%                                      =banwidth for 'flat' kernel;
%              flagUseKDTree: 0/1/2
%                             specifies which kdtree to use, if any:
%                             0 - no kdtree
%                             1 - matlab builtin kdtree
%                             2 - external kdtree (Andrea)
%                             Default: 1 (but recommended to use 2 if available, faster than matlab built-in but has problems in r2013b under windows)
%                            
%                  flagDebug: true/false
%                             specifies whether or not to run in debug mode.
%                             In debug mode, a bunch of stuff is/will-be printed 
%                             and plotted for debugging purposes.
%                            
%     Output Arguments:
%                      
%                clusterInfo: A structure array with one entry per cluster found
%                             The structure contains the following elements:
%                                 numPoints - number of points in the cluster
%                                 ptIdData - Ids (row indices of ptData) of 
%                                            the points belonging to the cluster
%                                 ptClusterCenter - coordinates of the cluster center
%                                 
%                                 
%          pointToClusterMap: vector of size equal to the number of points 
%                             in the input point set.
%                             pointToClusterMap(i) - cluster label of point i.   
%                
%         
% Dependencies:
% 
%       This code uses the kdtree implementation of Andrea Tagliasacchi 
%       (http://www.mathworks.com/matlabcentral/fileexchange/21512-kd-tree-for-matlab)
% 
% Author: Deepak Roy Chittajallu (modified from MeanShiftClustering by
% Hunter Elliott)
% 
% References:
% 
% 1) Comaniciu, D. and P. Meer (2002). 
%    "Mean shift: a robust approach toward feature space analysis.", 
%    IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%    24(5): 603-619.
% 

    p = inputParser;
    p.addRequired( 'ptData', @(x) (isnumeric(x) && ismatrix(x)) );
    p.addRequired( 'bandwidth', @(x) (size(x,ndims(x)) == size(ptData,1) || numel(x) == size(ptData,1))  );
    p.addParamValue( 'kernel', 'gaussian', @(x) ( (ischar(x) && ismember(x, {'gaussian', 'flat'})) || (isa(x,'function_handle') && nargin(x) >= 3 && abs(nargout(x)) >= 1)) );
    p.addParamValue( 'method', 'standard', @(x) ( (ischar(x) && ismember(x, {'standard', 'optimized'})) ) );
    p.addParamValue( 'maxIterations', 500, @(x) (isscalar(x)) );
    p.addParamValue('kernelSupport',[],@(x)(numel(x) == size(ptData,1)));
    p.addParamValue( 'minClusterDistance', .05,  @(x) (numel(x) == 1));
    p.addParamValue( 'flagUseKDTree', 1, @(x) (isscalar(x) && ismember(x,0:2)) );
    p.addParamValue( 'flagDebug', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( ptData, bandwidth, varargin{:} );
    
    kernelfunc = p.Results.kernel;
    meanShiftMethod = p.Results.method;
    maxIterations = p.Results.maxIterations;
    minClusterDistance = p.Results.minClusterDistance;      
    flagDebug =  p.Results.flagDebug;
    flagUseKDTree = p.Results.flagUseKDTree;
    
    if ischar(p.Results.kernel)
       
        switch p.Results.kernel        
            case 'gaussian'                
                kernelfunc = @update_mean_gaussian_kernel;
            case 'flat'                             
                error('This function only supports gaussian kernels!')
                %kernelfunc = @update_mean_flat_kernel;                
        end
                
    else
          kernelfunc = p.Results.kernel;          
    end
    
    %If scalar bandwidth was input, convert to diagonal matrix
    d = size(ptData,2);%Dimensionality
    n = size(ptData,1);
    if numel(bandwidth) == n
        %TEMP - way to vectorize this?
        tmp = bandwidth;
        bandwidth = zeros(d,d,n);
        for j = 1:d
            bandwidth(j,j,:) = tmp;
        end       
    end
        
    %Pre-compute bandwidth determinants and inverses to save time while
    %iterating
    %TEMP - way to vectorize this???
    bandDet = zeros(n,1);
    bandInv = zeros(d,d,n);    
    for i = 1:n
        bandDet(i)= det(bandwidth(:,:,i));
        bandInv(:,:,i) = inv(bandwidth(:,:,i));                
    end

    if isempty(p.Results.kernelSupport)        
        kernelSupport = 2*sqrt(bandwidth);        
    else
        kernelSupport = p.Results.kernelSupport;
    end

    switch meanShiftMethod
        
        case 'standard'
            
            [clusterInfo, pointToClusterMap,pointTraj] = StandardMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport,bandDet,bandInv);
            
        case 'optimized'
            error('This function doesn''t yet support the optimized algorithm!')
%             [clusterInfo, pointToClusterMap] = OptimizedMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport);            
%             pointTraj = {};
    end
    
end

function [clusterInfo, pointToClusterMap] = OptimizedMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport)

    threshConvergence = 1e-3 * bandwidth;    
    [numDataPoints, numDataDims] = size( ptData );           
    
    %% Run mean-shift for each data point and find its mode
    fprintf( 1, '\nMean-shift clustering on a dataset of %d points ...\n', numDataPoints );
    pointClusterVotes = zeros( numDataPoints, 1 );
    flagPointVisited = false( numDataPoints, 1 );
    numPointsProcessed = 0;
    ptIdRandomized = randperm(numDataPoints);
    clusterInfo = [];
    
        % build kd-tree over all points if requested  -- O(n log(n) log(n)) 
        if flagUseKDTree
            
            if flagDebug
                fprintf( 1, '\n\tBuilding KD-Tree ...\n' );        
            end                 
            if flagUseKDTree == 1            
                kdtree_points = KDTreeSearcher( ptData );  
            else
                kdtree_points = kdtree_build(ptData);
            end
        end
        
        % for each data point climb the non-parametric density function to find its mode
        meanIterationsElapsed = 0;
        if flagDebug
            fprintf( 1, '\n\tFinding the mode for each data point ...\n' );     
        end       
        
        for i = 1:numDataPoints

            curPtId = ptIdRandomized(i);
            
            if flagPointVisited(curPtId)
                continue;
            end
            
            numIterationsElapsed = 0;        
            ptOldMean = ptData(curPtId,:);

            flagPointOnPath = false( numDataPoints, 1 );
            flagPointOnPath( curPtId ) = true;
            
            while true            

                % get nearest points
                if flagUseKDTree == 1
                    [ptIdNearest] = rangesearch(kdtree_points, ptOldMean, max(max(kernelSupport(:,:,i))));
                    ptIdNearest = ptIdNearest{1};  
                elseif flagUseKDTree == 2
                    ptIdNearest = kdtree_ball_query(kdtree_points, ptOldMean, max(max(kernelSupport(:,:,i))));                       
                else
                    ptIdNearest = exhaustive_ball_query( ptData, ptOldMean, max(max(kernelSupport(:,:,i))));                       
                end
                
                ptNearest = ptData( ptIdNearest, : );                            
                flagPointOnPath( ptIdNearest ) = true;                                
                
                % call kernel to shift mean
                [ptNewMean] = kernelfunc( ptOldMean, ptNearest, bandwidth(ptIdNearest) ); 
                numIterationsElapsed = numIterationsElapsed + 1;

                % check for convergence
                if numIterationsElapsed >= maxIterations || norm(ptNewMean - ptOldMean) < threshConvergence(i)                
                    meanIterationsElapsed = meanIterationsElapsed + numIterationsElapsed;
                    break;
                else
                    ptOldMean = ptNewMean;
                end            

            end               
            
            ptClusterCenter = ptNewMean;
            
            % make a note of all the points that were visited on the path to this mode 
            flagPointVisited = flagPointVisited | flagPointOnPath;            
            
            % check if a close cluster is present already
            blnCloseClusterFound = false;
            closestClusterDist = [];
            closestClusterId = [];
            for cid = 1:numel( clusterInfo )
            
                curClusterDist = norm( ptClusterCenter - clusterInfo(cid).ptClusterCenter );
                if curClusterDist < minClusterDistance(i) 
                    
                    if ~blnCloseClusterFound || (blnCloseClusterFound && curClusterDist < closestClusterDist)                      
                        blnCloseClusterFound = true;
                        closestClusterId = cid;
                        closestClusterDist = curClusterDist;
                    end
                    
                end
                
            end
               
            if ~blnCloseClusterFound
                clusterInfo(end+1).ptClusterCenter = ptClusterCenter;
                pointClusterVotes( flagPointOnPath, numel(clusterInfo) ) = 1;                
            else
                clusterInfo(closestClusterId).ptClusterCenter = 0.5 * (ptClusterCenter + clusterInfo(closestClusterId).ptClusterCenter);
                pointClusterVotes( flagPointOnPath, closestClusterId ) = pointClusterVotes( flagPointOnPath, closestClusterId ) + 1;
            end            
            
            numPointsProcessed = numPointsProcessed + 1;
            
        end        
    
        meanIterationsElapsed = meanIterationsElapsed / numPointsProcessed;
        
        fprintf( 1, '\n\t%d clusters were found ...\n',  numel( clusterInfo ) ); 
        fprintf( 1, '\n\tMode seeking for each data took an average of %d iterations ...\n', round( meanIterationsElapsed ) ); 

        % assign each point to maximum voting cluster
        [ maxvote, pointToClusterMap] = max( pointClusterVotes, [], 2 );
        
    %% post process clusters
    for i = 1:numel( clusterInfo )
       clusterInfo(i).ptIdData = find( pointToClusterMap == i );
       clusterInfo(i).numPoints = numel( clusterInfo(i).ptIdData );
    end
end

function [clusterInfo, pointToClusterMap, pointTraj] = StandardMeanShift( ptData, bandwidth, kernelfunc, maxIterations, minClusterDistance, flagUseKDTree, flagDebug, kernelSupport,bandDet,bandI)

    threshConvergence = 1e-2 * minClusterDistance;
    [numDataPoints, numDataDims] = size( ptData );           
        
    pointTraj = cell(numDataPoints,1);            
    
    
    
    %% Run mean-shift for each data point and find its mode
    fprintf( 1, '\nMean-shift clustering on a dataset of %d points ...\n', numDataPoints );
    pointToClusterMap = zeros( numDataPoints, 1 );
    clusterInfo = [];
    
        % build kd-tree over all points if requested  -- O(n log(n) log(n)) 
        if flagUseKDTree
            
            if flagDebug
                fprintf( 1, '\n\tBuilding KD-Tree ...\n' );        
            end        
            if flagUseKDTree == 1            
                kdtree_points = KDTreeSearcher( ptData );  
            elseif flagUseKDTree == 2
                kdtree_points = kdtree_build(ptData);
            end
        end        
        
        % for each data point climb the non-parametric density function to find its mode
        meanIterationsElapsed = 0;
        if flagDebug
            fprintf( 1, '\n\tFinding the mode for each data point ...\n' );        
        end       
        
        for i = 1:numDataPoints

            numIterationsElapsed = 0;        
            ptOldMean = ptData(i,:);
            
            pointTraj{i} = nan(maxIterations,numDataDims);
            pointTraj{i}(1,:) = ptOldMean;

            while true            

                % get nearest points
                if flagUseKDTree == 1
                    [ptIdNearest] = rangesearch(kdtree_points, ptOldMean, max(max(kernelSupport(:,:,i))));
                    ptIdNearest = ptIdNearest{1};  
                elseif flagUseKDTree == 2
                    ptIdNearest = kdtree_ball_query(kdtree_points, ptOldMean, max(max(kernelSupport(:,:,i))));                       
                else
                    ptIdNearest = exhaustive_ball_query( ptData, ptOldMean, max(max(kernelSupport(:,:,i))));                       
                end
                ptNearest = ptData( ptIdNearest, : );
                
                % call kernel to shift mean
                [ptNewMean] = kernelfunc( ptOldMean, ptNearest, bandwidth(:,:,ptIdNearest),bandDet(ptIdNearest),bandI(:,:,ptIdNearest)); 
                numIterationsElapsed = numIterationsElapsed + 1;
                
                pointTraj{i}(numIterationsElapsed+1,:) = ptNewMean;

                % check for convergence using mahalanobis distance  
                dispVec = ptNewMean - ptOldMean;
                if numIterationsElapsed >= maxIterations || chi2cdf(dispVec * bandwidth(:,:,i) * dispVec',numDataDims) < threshConvergence
                    meanIterationsElapsed = meanIterationsElapsed + numIterationsElapsed;
                    break;
                else
                    ptOldMean = ptNewMean;
                end            

            end               

            ptClusterCenter = ptNewMean;     
            pointTraj{i} = pointTraj{i}(1:numIterationsElapsed+1,:);
            
            % get point density around the cluster center            
            if flagUseKDTree == 1
                [ptIdNearest] = rangesearch(kdtree_points, ptOldMean, max(max(kernelSupport(:,:,i))));
                ptIdNearest = ptIdNearest{1};  
            elseif flagUseKDTree == 2
                ptIdNearest = kdtree_ball_query(kdtree_points, ptOldMean, max(max(kernelSupport(:,:,i))));                       
            else
                ptIdNearest = exhaustive_ball_query( ptData, ptOldMean, max(max(kernelSupport(:,:,i))));                       
            end
            curClusterCenterDensity = numel( ptIdNearest );
            
            % check if cluster is present already
            blnCloseClusterFound = false;
            if ~isempty(clusterInfo)
                dX = bsxfun(@minus,ptClusterCenter,vertcat(clusterInfo(:).ptClusterCenter));%Separation vector for each cluster
                nC = numel(clusterInfo);
                Dm = zeros(nC,1);
                for j= 1:nC
                    Dm(j) = dX(j,:) * clusterInfo(j).clusterMeanHi * dX(j,:)';%Mahalanobis distance to each cluster center
                end
                                
                [minD,cid] = min(Dm);

                %Test distance based on elliptical estimate of percentile
                if chi2cdf(minD^2,numel(ptClusterCenter)) < minClusterDistance
                    blnCloseClusterFound = true;
                    if curClusterCenterDensity > clusterInfo(cid).clusterCenterDensity
                        clusterInfo(cid).ptClusterCenter = ptClusterCenter;
                        clusterInfo(cid).clusterCenterDensity = curClusterCenterDensity;                                                
                        [~,Hhi] = kernelfunc(ptClusterCenter,ptData(ptIdNearest,:), bandwidth(:,:,ptIdNearest),bandDet(ptIdNearest),bandI(:,:,ptIdNearest)); %Get harmonic mean of bandwidths for central points. TEMP - use all points??
                        clusterInfo(end).clusterMeanHi = Hhi;                                                
                    end
                    pointToClusterMap(i) = cid;

                end
            end
                        
            if ~blnCloseClusterFound
                clusterInfo(end+1).ptClusterCenter = ptClusterCenter;
                clusterInfo(end).clusterCenterDensity = curClusterCenterDensity;
                [~,Hhi] = kernelfunc(ptClusterCenter,ptData(ptIdNearest,:), bandwidth(:,:,ptIdNearest),bandDet(ptIdNearest),bandI(:,:,ptIdNearest)); %Get harmonic mean of bandwidths for central points. TEMP - use all points??
                clusterInfo(end).clusterMeanHi = Hhi;
                pointToClusterMap(i) = numel( clusterInfo );
            end            
            
        end        
    
        meanIterationsElapsed = meanIterationsElapsed / numDataPoints;
        
        fprintf( 1, '\n\t%d clusters were found ...\n',  numel( clusterInfo ) ); 
        fprintf( 1, '\n\tMode seeking for each data took an average of %d iterations ...\n', round( meanIterationsElapsed ) ); 
        
    %% post process clusters
    for i = 1:numel( clusterInfo )
       clusterInfo(i).ptIdData = find( pointToClusterMap == i );
       clusterInfo(i).numPoints = numel( clusterInfo(i).ptIdData );
    end

end


function [ptNewMean,Hhi] = update_mean_gaussian_kernel( ptOldMean, ptNearest, H,detH,Hi)

    %TEMP - vectorize ALL of this????
    dX = bsxfun(@minus,ptOldMean,ptNearest);
    n = size(ptNearest,1);
    d = size(ptNearest,2);
    Dm = zeros(n,1);    
    HidX = zeros(d,n);    
    for i = 1:n                
        HidX(:,i) = Hi(:,:,i) * dX(i,:)';%Store so we can reuse
        Dm(i) = dX(i,:) * HidX(:,i);%Mahalanobis distance to each point        
    end
    
    w = 1 ./ sqrt(detH) .* exp(-.5 * Dm);%Weight vector
    w = w ./ sum(w);
        
    Hhi = sum(bsxfun(@times,Hi,permute(w,[3 2 1])),d+1);%Inverse weighted harmonic mean of bandwidth matrices
            
    ptNewMean = ptOldMean - (Hhi \ sum(bsxfun(@times,HidX',w),1)')';%Mean-shifted position        
    
end


function [ptNewMean] = update_mean_flat_kernel( ptOldMean, ptNearest, bandwidth )

    ptNewMean = mean( ptNearest, 1 );
    
end

function [ ptIdNearest,dist] = exhaustive_ball_query( ptData, ptPoint, bandwidth )

    numDataPoints = size( ptData, 1 );
    sqDistances = sum( ( ptData - repmat( ptPoint, numDataPoints, 1 ) ).^2, 2 );
    ptIdNearest = find( sqDistances < bandwidth^2 );      
    dist = sqrt(sqDistances(ptIdNearest));
    
end