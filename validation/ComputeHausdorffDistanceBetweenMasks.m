function [ dist ] = ComputeHausdorffDistanceBetweenMasks( mask1 , mask2 , ImageSpacing )

    if ~exist( 'ImageSpacing' , 'var' )
        ImageSpacing = ones( ndims( mask1 ) , 1 );
    end
    
	mask1_perim = bwperim( mask1 > 0 , 8 );
	mask2_perim = bwperim( mask2 > 0 , 8 );
	
	mask1_dist = bwdistsc(mask1_perim,ImageSpacing);
	mask2_dist = bwdistsc(mask1_perim,ImageSpacing);
	
	dist = max( [ max( mask1_dist(mask2_perim > 0) ) , max( mask2_dist(mask1_perim > 0) ) ] );
	
%     bnd_pts_mask1 = GetPointsOnMaskBoundary( mask1 );
%     bnd_pts_mask2 = GetPointsOnMaskBoundary( mask2 );
% 
%     dist = hausdorff( bnd_pts_mask1, bnd_pts_mask2 , ImageSpacing );

end

function [ sub_ind , lin_ind ] = GetPointsOnMaskBoundary( mask )
   
    lin_ind = find( bwperim( mask > 0 , 8 ) );    
    sub_ind = cell( 1 , ndims( mask ) );    
    [ sub_ind{:} ] = ind2sub( size( mask ) , lin_ind );     
    sub_ind = [ sub_ind{:} ];
    
end