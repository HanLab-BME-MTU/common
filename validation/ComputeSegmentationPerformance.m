function [stats] = ComputeSegmentationPerformance( imSegMask , imTruthMask , ImageSpacing )

    if ~exist( 'ImageSpacing' , 'var' )

        ImageSpacing = ones( 1 , ndims(imSegMask) );

    end

    stats = ComputeClassificationPerformance( imSegMask(:) > 0 , imTruthMask(:) > 0 ); 

        % Additional stats
        stats.ImageSpacing = ImageSpacing;
        stats.pixel_area = prod( ImageSpacing(1:2) );

        stats.overlap_area = numel( find( imSegMask > 0 & imTruthMask > 0 ) );
        stats.overlap_area_physp = stats.overlap_area * stats.pixel_area;

        stats.seg_area = numel( find( imSegMask > 0 ) );
        stats.seg_area_physp = stats.seg_area * stats.pixel_area;

        stats.truth_area = numel( find( imTruthMask > 0 ) );
        stats.truth_area_physp = stats.truth_area * stats.pixel_area;           

        stats.DiceCoefficient = ( 2 * stats.overlap_area ) / ( stats.seg_area + stats.truth_area );     
        
        stats.HausdorffDistance = ComputeHausdorffDistanceBetweenMasks( imSegMask , imTruthMask , ImageSpacing );
        
end