function [ WAcc ] = ComputeWeightedAccuracy( stats )

    WAcc = zeros( numel( stats ) , 1 );
    
    for i = 1:numel( stats )
        
        PositiveRecall = stats(i).TP / ( stats(i).TP + stats(i).FN );
        NegativeRecall = stats(i).TN / ( stats(i).TN + stats(i).FP );

        WAcc(i) = mean( [ PositiveRecall NegativeRecall ] );
        
    end
    
end