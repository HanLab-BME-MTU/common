function [ PositiveRecall ] = ComputePositiveRecall( stats )

    PositiveRecall = zeros( numel( stats ) , 1 );
    
    for i = 1:numel( stats )
        
        PositiveRecall(i) = stats(i).TP / ( stats(i).TP + stats(i).FN );

    end

end