function [ boxsum ] = ComputeImageBoxIntegral( imIntegral, boxOrigin, boxSize )
    
    imdims = ndims(imIntegral);    

    switch imdims
        
        case 2
        
            boxsum = imIntegral( boxOrigin(1), boxOrigin(2) ) ...
                     - imIntegral( boxOrigin(1), boxOrigin(2) + boxSize(2) ) ... 
                     - imIntegral( boxOrigin(1) + boxSize(1), boxOrigin(2) ) ...
                     + imIntegral( boxOrigin(1) + boxSize(1), boxOrigin(2) + boxSize(2) );
                 
        case 3
            
            boxsum = - imIntegral( boxOrigin(1), boxOrigin(2), boxOrigin(3) ) ...
                     + imIntegral( boxOrigin(1), boxOrigin(2) , boxOrigin(3) + boxSize(3) ) ...
                     + imIntegral( boxOrigin(1), boxOrigin(2) + boxSize(2), boxOrigin(3) ) ...
                     - imIntegral( boxOrigin(1), boxOrigin(2) + boxSize(2), boxOrigin(3) + boxSize(3) ) ...
                     + imIntegral( boxOrigin(1) + boxSize(1), boxOrigin(2), boxOrigin(3) ) ...
                     - imIntegral( boxOrigin(1) + boxSize(1), boxOrigin(2), boxOrigin(3) + boxSize(3) ) ...
                     - imIntegral( boxOrigin(1) + boxSize(1), boxOrigin(2) + boxSize(2), boxOrigin(3) ) ...
                     + imIntegral( boxOrigin(1) + boxSize(1), boxOrigin(2) + boxSize(2), boxOrigin(3) + boxSize(3) );                
        
        otherwise
            
            imsize = size( imIntegral );
            numCorners = 2^imdims;
            boxsum = 0;   
    
            for i = 0:numCorners-1
                b = dec2binvec(i,ndims(imIntegral));                
                cornerPixSubind = num2cell( boxOrigin(:) + b' .* boxSize(:) );       
                cornerPixInd = sub2ind( imsize, cornerPixSubind{:} );
                boxsum = boxsum + (-1)^(imdims-sum(b)) * imIntegral( cornerPixInd );
            end
        
    end
    
end

