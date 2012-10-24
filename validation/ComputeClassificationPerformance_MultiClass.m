function [ stats ] = ComputeClassificationPerformance_MultiClass( PredictedLabels , ActualLabels , LabelList )
% 
%                             Predicted
%                     class1  class2  class3 ... ... 
%         class1
%         class2
% Actual  class3
%         ...
%         ...
% 

    if ~exist( 'LabelList' , 'var' )
        
        LabelList = unique( cat( 1 , PredictedLabels , ActualLabels ) );
        
    end
    
    stats.LabelList = LabelList;
    stats.PredictedLabels = PredictedLabels;
    stats.ActualLabels = ActualLabels;
    
    numLabels = numel( LabelList );
    
    % compute multi-class confusion matrix
    stats.cMatrix = zeros( numLabels , numLabels );
    
    for i = 1:numLabels
        
        for j = 1:numLabels
        
            stats.cMatrix( i , j ) = numel( find( ActualLabels == LabelList(i) & PredictedLabels == LabelList(j) ) );           
            
        end
        
    end        

    stats.PredictionAccuracy = sum( diag( stats.cMatrix ) ) * 100 / numel( PredictedLabels );   

    stats.cMatrix_Display = zeros( size( stats.cMatrix ) + 1 );    
    stats.cMatrix_Display( 2:end , 2:end ) = stats.cMatrix;
    stats.cMatrix_Display( 2:end , 1 ) = stats.LabelList;
    stats.cMatrix_Display( 1 , 2:end ) = stats.LabelList;
    stats.cMatrix_Display( 1 , 1) = -1;
    
    % compute per-label stats
    for i = 1:numLabels
        
        curPredictedLabels = PredictedLabels;
        curPredictedLabels( PredictedLabels == LabelList(i) ) = 1;
        curPredictedLabels( PredictedLabels ~= LabelList(i) ) = 0;

        curActualLabels = ActualLabels;
        curActualLabels( ActualLabels == LabelList(i) ) = 1;        
        curActualLabels( ActualLabels ~= LabelList(i) ) = 0;

        stats.PerClassStats(i) = ComputeClassificationPerformance( curPredictedLabels > 0 , curActualLabels > 0 );
        
    end
    
end