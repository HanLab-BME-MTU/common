% Francois Aguet, 10/24/2012

function CC = computeSegmentProperties(CC, img, theta)

dims = size(img);
ny = dims(1);
[x,y] = meshgrid(1:dims(2),1:ny);

getNH = @(i) [i-ny-1 i-ny i-ny+1 i-1 i+1 i+ny-1 i+ny i+ny+1];

ns = CC.NumObjects;
CC.rawAngle = cell(1,ns);
CC.rval = cell(1,ns);
CC.lval = cell(1,ns);
for i = 1:ns
    
    np = numel(CC.PixelIdxList{i});
    unorderedIdx = CC.PixelIdxList{i};
    orderedIdx = zeros(np,1);
    
    % assign 1st endpoint (choice is arbitrary)
    orderedIdx(1) = CC.endpointIdx{i}(1);
    unorderedIdx(unorderedIdx==orderedIdx(1)) = [];
    
    for k = 2:np
        % local neighborhood indexes
        hoodIdx = getNH(orderedIdx(k-1));
        nextIdx = intersect(unorderedIdx, hoodIdx);
        if isempty(nextIdx)
            % get endpoints of unordered index
            umask = zeros(dims);
            umask(unorderedIdx) = 1;
            %umask = double(bwmorph(umask, 'thin'));
            nn = (imfilter(umask, ones(3), 'same')-1) .* umask;
            endpointIdx = find(nn<2 & umask==1);
            
            [yi, xi] = ind2sub(dims, orderedIdx(orderedIdx~=0));
            X = [xi yi];
            [yi, xi] = ind2sub(dims, endpointIdx);
            x0 = [xi yi];
            [~,dist] = KDTreeClosestPoint(X, x0);
            
            nextIdx = endpointIdx(find(dist==min(dist),1,'first')); % start over at closest endpoint
        end
        unorderedIdx(unorderedIdx==nextIdx(1)) = [];
        orderedIdx(k) = nextIdx(1);
    end
    CC.PixelIdxList{i} = orderedIdx;
    CC.rawAngle{i} = theta(orderedIdx);
    cost = cos(CC.rawAngle{i});
    sint = sin(CC.rawAngle{i});
    [yi, xi] = ind2sub(dims, CC.PixelIdxList{i});
    X = [xi+cost xi+2*cost];
    Y = [yi+sint yi+2*sint];
    CC.rval{i} = interp2(x, y, img, X, Y);
    X = [xi-cost xi-2*cost];
    Y = [yi-sint yi-2*sint];
    CC.lval{i} = interp2(x, y, img, X, Y);
    
    % running average - POTENTIAL BUG: angle wrapping in (-pi/2..pi/2]
    %if np>w
    %    CC.smoothAngle{i} = conv([CC.rawAngle{i}(b:-1:2) CC.rawAngle{i} CC.rawAngle{i}(end-1:-1:end-b+1)], ones(1,w)/w, 'valid');
    %else
    %    CC.smoothAngle{i} = CC.rawAngle{i};
    %end
end

