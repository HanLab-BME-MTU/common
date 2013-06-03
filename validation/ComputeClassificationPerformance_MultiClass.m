function [ stats ] = ComputeClassificationPerformance_MultiClass( inPredictedLabels , inActualLabels , inLabelList )
% 
%                             Predicted
%                     class1  class2  class3 ... ... 
%         class1
%         class2
% Actual  class3
%         ...
%         ...
% 

    if ~exist( 'inLabelList' , 'var' )        
        inLabelList = unique( cat( 1 , inPredictedLabels , inActualLabels ) );
    end
    
    stats.LabelList = inLabelList;
    stats.PredictedLabels = inPredictedLabels;
    stats.ActualLabels = inActualLabels;
    
    if iscell(inPredictedLabels)
        
        PredictedLabels = zeros( size(inPredictedLabels) );
        ActualLabels = zeros( size(inPredictedLabels) );
        
        for i = 1:numel(PredictedLabels)
            curPredictedLabel = find( strcmp( inPredictedLabels{i}, inLabelList ) );
            curActualLabel = find( strcmp( inActualLabels{i}, inLabelList ) );
            
            if isempty(curPredictedLabel) || isempty(curActualLabel)
               error( 'error: found a label which is not in the label list' ); 
            end
            
            PredictedLabels(i) = curPredictedLabel;
            ActualLabels(i) = curActualLabel;
        end
        
        LabelList = 1:numel(inLabelList);
        
    else        
        LabelList = inLabelList;
        PredictedLabels = inPredictedLabels;
        ActualLabels = inActualLabels;
    end
    
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
    stats.cMatrix_Display( 2:end , 1 ) = LabelList;
    stats.cMatrix_Display( 1 , 2:end ) = LabelList;
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
    
    stats.GMean = geomean( [stats.PerClassStats.Recall] );
    stats.avgFMeasure = mean( [stats.PerClassStats.FMeasure] );
    
end