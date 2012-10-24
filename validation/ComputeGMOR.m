function [ GMOR ] = ComputeGMOR( stats )

    GMOR = zeros( numel( stats ) , 1 );
    
    for i = 1:numel( stats )
        
        PositiveRecall = stats(i).TP / ( stats(i).TP + stats(i).FN );
        NegativeRecall = stats(i).TN / ( stats(i).TN + stats(i).FP );

        GMOR(i) = geomean( [ PositiveRecall NegativeRecall ] );
        
    end
    
end