function [imLabelRGB, varargout] = label2rgbND( imLabel, labelColorMap )

    if ~exist( 'labelColorMap', 'var' )

        % assign colors randomly
        numLabels = double(max(imLabel(:)));
        imLabel( imLabel == 0 ) = numLabels + 1;
        cmap = jet(numLabels);
        stream = RandStream('swb2712','seed',0);
        index = randperm(stream,numLabels);
        cmap = [ cmap(index,:); 0 0 0 ];            
        
    else
        assert( all(size(labelColorMap, 2) == 3) );
        assert( all(size(labelColorMap, 1) >= max(imLabel(:))+1) );                
        
        numLabels = size(labelColorMap, 1) - 1; 
        imLabel( imLabel == 0 ) = numLabels + 1;
        
        cmap = labelColorMap;
    end
    
    imLabelRGB = [];
    for i = 1:3
        imLabelRGB = cat( ndims(imLabel) + 1, imLabelRGB, reshape( cmap( imLabel(:), i ), size( imLabel ) ) );
    end
    
    if nargout > 1
        varargout{1} = cmap;
    end
    
end