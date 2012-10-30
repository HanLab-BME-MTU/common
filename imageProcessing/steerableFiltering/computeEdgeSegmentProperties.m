% order the pixels of each segment from one endpoint to the other

% Francois Aguet, 10/24/2012

function CC = computeEdgeSegmentProperties(CC, img, theta, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('InterpDistance', [1 2]); % [1 2 9 10]
ip.parse(varargin{:});
R = ip.Results.InterpDistance;
nr = numel(R);

[ny,nx] = size(img);
pixelIdxList = vertcat(CC.PixelIdxList{:});

segmentMask = zeros(ny,nx);
segmentMask(pixelIdxList) = 1;
labels = double(labelmatrix(CC));

% Identify endpoints
nn = (imfilter(segmentMask, ones(3), 'same')-1) .* segmentMask;
endpointIdxList = pixelIdxList(nn(pixelIdxList)==1);
% Each segment must have two endpoints
endpointIdxList = [endpointIdxList(1:2:end) endpointIdxList(2:2:end)];
tmp = NaN(CC.NumObjects,2);
tmp(labels(endpointIdxList(:,1)),:) = endpointIdxList;
endpointIdxList = tmp;

% Run a geodesic distance transform to calculate the pixel order
tmp = endpointIdxList(:,1);
D = bwdistgeodesic(logical(segmentMask), tmp(~isnan(tmp)));
D(isinf(D)) = 0;
CC.PixelOrder = mat2cell(D(pixelIdxList)+1, CC.NumPixels, 1);
CC.EndpointIdx = mat2cell(endpointIdxList, ones(size(endpointIdxList,1),1), 2);
for i = 1:CC.NumObjects
    CC.PixelIdxList{i} = CC.PixelIdxList{i}(CC.PixelOrder{i});
end

% update
pixelIdxList = vertcat(CC.PixelIdxList{:});

% compute intensity on side of all segments
angleVect = theta(pixelIdxList);
cost = cos(angleVect);
sint = sin(angleVect);
[x,y] = meshgrid(1:nx,1:ny);
[yi, xi] = ind2sub([ny nx], pixelIdxList);

% 'positive' or 'right' side
% X = [xi+R(1)*cost xi+R(2)*cost];
X = repmat(xi, [1 nr])+cost*R;
% Y = [yi+R(1)*sint yi+R(2)*sint];
Y = repmat(yi, [1 nr])+sint*R;
interp2(x, y, img, X, Y);
CC.rval = mat2cell(interp2(x, y, img, X, Y), CC.NumPixels, nr);

% 'negative' or 'left' side
% X = [xi-cost xi-2*cost];
% Y = [yi-sint yi-2*sint];
X = repmat(xi, [1 nr])-cost*R;
Y = repmat(yi, [1 nr])-sint*R;
CC.lval = mat2cell(interp2(x, y, img, X, Y), CC.NumPixels, nr);




% getNH = @(i) [i-ny-1 i-ny i-ny+1 i-1 i+1 i+ny-1 i+ny i+ny+1];
% 
% ns = CC.NumObjects;
% CC.rawAngle = cell(1,ns);
% CC.rval = cell(1,ns);
% CC.lval = cell(1,ns);
% for i = 1:ns
%     
%     np = numel(CC.PixelIdxList{i});
%     unorderedIdx = CC.PixelIdxList{i};
%     orderedIdx = zeros(np,1);
%     
%     % assign 1st endpoint (choice is arbitrary)
%     orderedIdx(1) = CC.EndpointIdx{i}(1);
%     unorderedIdx(unorderedIdx==orderedIdx(1)) = [];
%     
%     for k = 2:np
%         % local neighborhood indexes
%         hoodIdx = getNH(orderedIdx(k-1));
%         nextIdx = intersect(unorderedIdx, hoodIdx);
%         if isempty(nextIdx)
%             % get endpoints of unordered index
%             umask = zeros(dims);
%             umask(unorderedIdx) = 1;
%             %umask = double(bwmorph(umask, 'thin'));
%             nn = (imfilter(umask, ones(3), 'same')-1) .* umask;
%             endpointIdx = find(nn<2 & umask==1);
% 
%             if ~isempty(endpointIdx)
%                 [yi, xi] = ind2sub(dims, orderedIdx(orderedIdx~=0));
%                 X = [xi yi];
%                 [yi, xi] = ind2sub(dims, endpointIdx);
%                 x0 = [xi yi];
%                 [~,dist] = KDTreeClosestPoint(X, x0);
%                 nextIdx = endpointIdx(find(dist==min(dist),1,'first')); % start over at closest endpoint
%             else
%                 nextIdx = unorderedIdx(1);
%             end
%         end
%         unorderedIdx(unorderedIdx==nextIdx(1)) = [];
%         orderedIdx(k) = nextIdx(1);
%     end
%     CC.PixelIdxList{i} = orderedIdx;
%     CC.rawAngle{i} = theta(orderedIdx);
%     cost = cos(CC.rawAngle{i});
%     sint = sin(CC.rawAngle{i});
%     [yi, xi] = ind2sub(dims, CC.PixelIdxList{i});
%     X = [xi+cost xi+2*cost];
%     Y = [yi+sint yi+2*sint];
%     CC.rval{i} = interp2(x, y, img, X, Y);
%     X = [xi-cost xi-2*cost];
%     Y = [yi-sint yi-2*sint];
%     CC.lval{i} = interp2(x, y, img, X, Y);
%     
%     % running average - POTENTIAL BUG: angle wrapping in (-pi/2..pi/2]
%     %if np>w
%     %    CC.smoothAngle{i} = conv([CC.rawAngle{i}(b:-1:2) CC.rawAngle{i} CC.rawAngle{i}(end-1:-1:end-b+1)], ones(1,w)/w, 'valid');
%     %else
%     %    CC.smoothAngle{i} = CC.rawAngle{i};
%     %end
% end

