function [outsideIdx] = maskVectors(dispMatX,dispMatY,bwstackImg)
Npoints = length(dispMatX);
outsideIdx = false(Npoints,1);

for ii=1:Npoints
    if bwstackImg(round(dispMatY(ii)),round(dispMatX(ii)))
        outsideIdx(ii) = true;
    end
end
