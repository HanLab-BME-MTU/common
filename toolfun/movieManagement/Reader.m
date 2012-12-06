classdef  Reader < handle
    % Concrete implementation of MovieObject for a single movie
    
    properties
        sizeX
        sizeY
        sizeZ
        sizeC
        sizeT
    end
    
    methods(Abstract)
        getSizeX(obj)
        getSizeY(obj)
        getSizeZ(obj)
        getSizeC(obj)
        getSizeT(obj)
        loadImage(obj)
    end
    
end