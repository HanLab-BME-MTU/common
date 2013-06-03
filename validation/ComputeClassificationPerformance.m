function [ stats ] = ComputeClassificationPerformance( PredictedLabels , ActualLabels )

    cMatrix  = getConfusionMatrix(PredictedLabels,ActualLabels);
    [totaccuracy tpr tnr fpf fnf] = getConfusionMatrixStats(cMatrix);

    TN = cMatrix(1,1);
    FN = cMatrix(2,1);
    FP = cMatrix(1,2);
    TP = cMatrix(2,2);                        
    
    stats.TP = TP;
    stats.FP = FP;
    stats.TN = TN;
    stats.FN = FN;

    stats.PredictionAccuracy = totaccuracy;
    stats.tpr = tpr;
    stats.tnr = tnr;
    stats.fpf = fpf;   
    stats.fnf = fnf;     

    Precision = TP / ( TP + FP );                        
    Recall = TP / ( TP + FN );             
    if isnan(Recall)
        Recall = 1;
    end
    
    if isnan(Precision)
        Precision = 1;
    end

    FMeasure = 2 * Precision * Recall / ( Precision + Recall );
    
    PositiveRecall = TP / ( TP + FN );
    NegativeRecall = TN / ( TN + FP );
    
    stats.Precision = Precision * 100;
    stats.Recall = Recall * 100;
    
    stats.PositiveRecall = PositiveRecall * 100;
    stats.NegativeRecall = NegativeRecall * 100;    
    
    stats.FMeasure = FMeasure * 100;              
    stats.GMOR = ComputeGMOR( stats ) * 100;
    stats.WAcc = ComputeWeightedAccuracy( stats ) * 100;
    
end