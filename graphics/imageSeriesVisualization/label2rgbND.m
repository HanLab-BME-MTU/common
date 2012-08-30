function [imLabelRGB, varargout] = label2rgbND( imLabel )

    numLabels = double(max(imLabel(:)));
    cmap = jet(numLabels);
    stream = RandStream('swb2712','seed',0);
    index = randperm(stream,numLabels);
    cmap = [ cmap(index,:); 0 0 0 ];        
    imLabel( imLabel == 0 ) = numLabels + 1;
    imLabelRGB = [];
    for i = 1:3
        imLabelRGB = cat( ndims(imLabel) + 1, imLabelRGB, reshape( cmap( imLabel(:), i ), size( imLabel ) ) );
    end
    
    if nargout > 1
        varargout{1} = cmap;
    end
    
end