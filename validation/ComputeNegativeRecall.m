function [ NegativeRecall ] = ComputeNegativeRecall( stats )

    NegativeRecall = zeros( numel( stats ) , 1 );
    
    for i = 1:numel( stats )
        
        NegativeRecall(i) = stats(i).TN / ( stats(i).TN + stats(i).FP );

    end

end