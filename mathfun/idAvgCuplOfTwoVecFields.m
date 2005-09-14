function S = idAvgCOfTwoVecFields(Va,Vh,P,varargin)
%idAvgCuplOfTwoVecFields: This function identifies the average coupling 
%                         between two vector fields.
%
% SYNOPSIS: 
%    S = idAvgCuplOfTwoVecFields(Va,Vh,P)
%    S = idAvgCuplOfTwoVecFields(Va,Vh,P,par1,value1, ...)
%
% INPUT:
%    Va : An n-by-2 matrix that stores data of the first vector field where n
%         is the number of vectors.
%    Vh : An n-by-2 matrix that stores data of the second vector field.
%     P : An n-by-2 matrix that store the spatial location of each vector.
%
% OPTIONAL Par/Value PAIRS:
%    'blockSize'      : The side length of the basic block used for
%                       subpopulation averaging and coherence test. 
%                       Default is 30 pixels.
%    'scaling'        : 'on' or 'off' (default). Specify whether the two
%                       vector fields are scaled before averging so that 
%                       different speed population can have the same weight 
%                       in calculating the coupling.
%    'smSpdThreshold' : Threshold for cutting off small speed that can not be
%                       trusted for statistic analysis. Default value is 0.
%    'cohrThreshold'  : Threshold for cutting off local (basic block) vector 
%                       population whose spatial coherence is too low to be
%                       trusted for statistic analysis. Default value is 0.05.
%    'pplThreshold'   : Population threshold for statistic analysis. Default
%                       value is 3;
%    'figH'           : The ID number of the figure window where the cluser 
%                       image of vector subpopulation for averaging is 
%                       displayed along with the two vector fields.
%    'image'          : An image (or image file) can also be provided for 
%                       covering the vector field.
%    'dispScale       : The scaling factor for drawing vectors.
%                       Default is [] meaning that no figure is shown.
%    'couplingVector' : For testing purpose, the unknown coupling vector can
%                       also be given (e.g from simulated data) so that all
%                       rotation will be based on this vector. Default is [].
%
% OUTPUT:
%    S : A structure that contains the following fields:
%        'alpha'     : the coupling coefficient.
%        'dAngl'     : the deviation angle between 'Vac' and 'Vhc'.
%        'aCohrS' : The cohrence score of 'Va'.
%        'hCohrS' : The cohrence score of 'Vh'.
%
% The Model:
%
%          Va = Vac       + Ra
%          Vh = alpha*Vhc + Rh
%
% where
%    Vac, Vhc : The coupling vectors between the two fields and |Vac| = |Vhc|. 
%               In an ideal coupling, Vac = Vhc. Otherwise, the angle between 
%               Vac and Vhc is defined as the deviation angle.
%    Ra       : Random vector (or noise) in the first field.
%    Rh       : Random vector (or noise) in the second field.
%    alpha    : Coupling coefficient. 
%
% Assumptions:
%    1. The means of 'Ra' and 'Rh' are zeros.
%
% Data:
%    meanVa : The mean of 'Va' after aligning all vectors by Vac.
%    meanVh : The mean of 'Vh' after aligning all vectors by Vac.
%    Da     : The norm of 'meanVa'.
%    Dh     : The norm of 'meanVh'.
%
% Small speed and coherence test for outlier exclusion:
%    A coherence score is calculated by the ratio between the length of the
%    average vector and the average speed (scaled if 'scaling' is 'on') for
%    each basic block.
%    A basic block is excluded if both the averge speed (not scaled) of the 
%    block population is below 'smSpdThreshold' and the coherence score is 
%    below 'cohrThreshold'.

%Default paramters.
bSize          = 30;
scaling        = 'on';
smSpdThreshold = 0;
cohrThreshold  = 0.05;
pplThreshold   = 3;
figH           = [];
img            = [];
Vc             = Va;
dispScale      = 15;

%Default return.
S.alpha  = NaN;
S.dAngl  = NaN;
S.aCohrS = NaN;
S.hCohrS = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parsing optional parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(nargin-3,2) ~= 0
   error('The optional parameter/value pair is not even.');
end

for k = 1:2:nargin-3
   switch varargin{k}
      case 'blockSize'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            bSize = varargin{k+1};
         end
      case 'scaling'
         scaling = varargin{k+1};
      case 'smSpdThreshold'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            smSpdThreshold = varargin{k+1};
         end
      case 'cohrThreshold'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            cohrThreshold = varargin{k+1};
         end
      case 'pplThreshold'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            pplThreshold = varargin{k+1};
         end
      case 'figH'
         figH = varargin{k+1};
      case 'image'
         img = varargin{k+1};
      case 'dispScale'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            dispScale = varargin{k+1};
         end
      case 'couplingVector'
         if ~isempty(varargin{k+1})
            Vc = varargin{k+1};
         end
      otherwise
         error(['Optional parameter '' varargin{k} '' is not recogonized.']);
   end
end

%Remove NaN from Va and Vh.
nanInd = find(isnan(Va(:,1)) | isnan(Va(:,2)) | ...
   isnan(Vh(:,1)) | isnan(Vh(:,2)) | isnan(Vc(:,1)) | isnan(Vc(:,2)));
Va(nanInd,:) = [];
Vh(nanInd,:) = [];
Vc(nanInd,:) = [];
P(nanInd,:)  = [];

if isempty(Va)
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create an image area that covers the vector field.
corLen = ceil(bSize/2);

cSpd = sqrt(sum(Vc.^2,2));
aSpd = sqrt(sum(Va.^2,2));
hSpd = sqrt(sum(Vh.^2,2));

%To avoid dividing by zero.
smNumber = min([1e-6 1e-6*median(aSpd) 1e-6*median(hSpd)]);
if smNumber == 0
   smNumber = 1e-6;
end

unitVa = NaN*ones(size(Va));
unitVh = NaN*ones(size(Vh));

nzInd = find(aSpd>=smNumber);
if ~isempty(nzInd)
   unitVa(nzInd,1) = Va(nzInd,1)./aSpd(nzInd);
   unitVa(nzInd,2) = Va(nzInd,2)./aSpd(nzInd);
end

nzInd = find(hSpd>=smNumber);
if ~isempty(nzInd)
   unitVh(nzInd,1) = Vh(nzInd,1)./hSpd(nzInd);
   unitVh(nzInd,2) = Vh(nzInd,2)./hSpd(nzInd);
end

%Note: Before rotation, rotVc (or rotVa or rotVh) means the scaled Vc if
% scaling is on. Otherwise, it is the same as Vc.
rotVc = Vc;
rotVa = Va;
rotVh = Vh;
if strcmp(scaling,'on')
   %'rotVc' is normalized 'Vc' used in calculating the direction for rotation.
   nzInd = find(cSpd>=smNumber);
   if ~isempty(nzInd)
      rotVc(nzInd,1) = Vc(nzInd,1)./cSpd(nzInd);
      rotVc(nzInd,2) = Vc(nzInd,2)./cSpd(nzInd);
   end

   %We also scale 'Va' and 'Vh' so that coupling of different speed population
   % are equally weighted. Note that we divide by 'aSpd+hSpd' to avoid
   % dividing by a very small number. Since the two channels are divided by
   % the same number, it will not affect the calculated average coupling.
   bgSpdInd = find(aSpd>smSpdThreshold);
   smSpdInd = find(aSpd<=smSpdThreshold);
   
   if ~isempty(bgSpdInd)
      rotVa(bgSpdInd,1) = Va(bgSpdInd,1)./(aSpd(bgSpdInd)+hSpd(bgSpdInd));
      rotVa(bgSpdInd,2) = Va(bgSpdInd,2)./(aSpd(bgSpdInd)+hSpd(bgSpdInd));
      rotVh(bgSpdInd,1) = Vh(bgSpdInd,1)./(aSpd(bgSpdInd)+hSpd(bgSpdInd));
      rotVh(bgSpdInd,2) = Vh(bgSpdInd,2)./(aSpd(bgSpdInd)+hSpd(bgSpdInd));
   end

   if ~isempty(smSpdInd)
      rotVa(smSpdInd,1) = Va(smSpdInd,1)./(smSpdThreshold+hSpd(smSpdInd));
      rotVa(smSpdInd,2) = Va(smSpdInd,2)./(smSpdThreshold+hSpd(smSpdInd));
      rotVh(smSpdInd,1) = Vh(smSpdInd,1)./(smSpdThreshold+hSpd(smSpdInd));
      rotVh(smSpdInd,2) = Vh(smSpdInd,2)./(smSpdThreshold+hSpd(smSpdInd));
   end
end

if isinf(corLen)
   %When corLen is Inf, there is only one block and we calculate the average
   % of the whole population. Do nothing here except for the speed, coherence
   % and population test.
   avgASpd = mean(aSpd);

   nInd = find(~isnan(unitVa(:,1)) & ~isnan(unitVa(:,2)));
   if ~isempty(nInd)
      S.aCohrS  = norm(mean(unitVa(nInd,:)))/length(nInd);
   else
      S.aCohrS  = NaN;
   end

   nInd = find(~isnan(unitVh(:,1)) & ~isnan(unitVh(:,2)));
   if ~isempty(nInd)
      S.hCohrS  = norm(mean(unitVh(nInd,:)))/length(nInd);
   else
      S.hCohrS  = NaN;
   end

   if (avgASpd <= smSpdThreshold & (isnan(S.aCohrS) | ...
      S.aCohrS <= cohrThreshold)) | ...
      length(aSpd) <= pplThreshold
      return;
   end
else
   if isempty(img)
      minX = min(P(:,1));
      maxX = max(P(:,1));
      minY = min(P(:,2));
      maxY = max(P(:,2));
      x = minX-corLen:2*corLen:maxX+3*corLen;
      y = minY-corLen:2*corLen:maxY+3*corLen;

      %Translate so that the image area starts with pixel (1,1).
      P(:,1) = P(:,1)-x(1)+1;
      P(:,2) = P(:,2)-y(1)+1;
      x = x - x(1) + 1;
      y = y - y(1) + 1;
      img = zeros(y(end),x(end));
   else
      if ischar(img)
         img = imread(img);
      end
      x = 1:2*corLen:size(img,2);
      y = 1:2*corLen:size(img,1);
      img(:) = 0;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Divide the image area into blocks and identify the vector population for
   %each block.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   nGridY  = length(y);
   nGridX  = length(x);
   nBlocks = (nGridX-1)*(nGridY-1);

   b.imgPixInd = cell(nBlocks,1);
   k = 0;
   for jx = 1:nGridX-1
      for jy = 1:nGridY-1
         k = k+1;
         [bX bY] = meshgrid(x(jx):x(jx+1)-1,y(jy):y(jy+1)-1);
         b.imgPixInd{k} = sub2ind(size(img),bY(:),bX(:));
      end
   end

   %Round P to integer pixel and get the linear index.
   rndP = round(P);
   indP = sub2ind(size(img),rndP(:,2),rndP(:,1));

   %anglVc = sign(Vc(:,2)).*acos(Vc(:,1)./cSpd);

   b.vecInd     = cell(nBlocks,1);
   b.vecSetID   = zeros(nBlocks,1);
   b.center     = zeros(nBlocks,2);
   b.angl       = NaN*ones(nBlocks,1);
   b.avgVc      = NaN*ones(nBlocks,2);
   b.avgCSpd    = NaN*ones(nBlocks,1);
   b.aCohrS     = NaN*ones(nBlocks,1);
   b.hCohrS     = NaN*ones(nBlocks,1);
   b.outlier    = [];
   b.inlier     = 1:nBlocks;
   for k = 1:nBlocks
      %Identify the index boundary of the block and find the index of vectors 
      % that are inside each block.
      [yI xI]     = ind2sub([nGridY-1,nGridX-1],k);
      b.center(k,:) = [(x(xI)+x(xI+1))/2 (y(yI)+y(yI+1))/2];
      b.vecInd{k} = find(rndP(:,2)>=y(yI) & rndP(:,2) < y(yI+1) & ...
         rndP(:,1)>=x(xI) & rndP(:,1)<x(xI+1));

      if ~isempty(b.vecInd{k})
         %The average angle, vector and speed for each block.
         nInd = b.vecInd{k}(find(~isnan(unitVa(b.vecInd{k},1)) & ...
            ~isnan(unitVa(b.vecInd{k},2))));
         if ~isempty(nInd)
            b.aCohrS(k) = norm(mean(unitVa(nInd,:)))/length(nInd);
         end

         nInd = b.vecInd{k}(find(~isnan(unitVh(b.vecInd{k},1)) & ...
            ~isnan(unitVh(b.vecInd{k},2))));
         if ~isempty(nInd)
            b.hCohrS(k) = norm(mean(unitVh(nInd,:)))/length(nInd);
         end
         b.avgASpd(k) = mean(sqrt(sum(unitVa(b.vecInd{k},:).^2,2)));
         
         %To avoid dividing by zero, we use 'smNumber' defined above.
         b.avgVc(k,:) = mean(Vc(b.vecInd{k},:));
         b.angl(k)    = sign(b.avgVc(k,2))* ...
                        acos(b.avgVc(k,1)/max(smNumber,norm(b.avgVc(k,:))));

         %Speed, coherence and population test.
         if (b.avgASpd(k) <= smSpdThreshold & ...
            (isnan(b.aCohrS(k)) | b.aCohrS(k) <= cohrThreshold)) | ...
            length(b.vecInd{k}) <= pplThreshold
            b.outlier = [b.outlier; k];
         end
      end
   end
   b.inlier([b.outlier; find(isnan(b.angl))]) = [];
   
   numberInd = find(~isnan(b.aCohrS));
   if ~isempty(numberInd)
      S.aCohrS = mean(b.aCohrS(numberInd));
   end
   numberInd = find(~isnan(b.hCohrS));
   if ~isempty(numberInd)
      S.hCohrS = mean(b.hCohrS(numberInd));
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Group blocks into a set of angle bins. We assign angle bins according the
   % natural angle distribution. This is a little bit tricky to code.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %Identify clusters of angle distribution.
   % First, get the population that are withing pi/6 (30 degree) from each
   % angle node in 'anglNodes' which is a vector of evenly distributed angle 
   % points between -pi and pi. 
   anglNodes = linspace(-pi,pi,73);
   anglDstr  = zeros(size(anglNodes));
   inlierAngl = b.angl(b.inlier);
   for k = 1:length(anglNodes)
      if anglNodes(k)-pi/6 < -pi
         anglDstr(k) = length(find((inlierAngl<=anglNodes(k)+pi/6 & ...
            inlierAngl>=-pi) | ...
            (inlierAngl<=pi & inlierAngl>=2*pi+anglNodes(k)-pi/6)));
      elseif anglNodes(k)+pi/6 > pi
         anglDstr(k) = length(find((inlierAngl>=anglNodes(k)-pi/6 & ...
            inlierAngl<=pi) | ...
            (inlierAngl>=-pi & inlierAngl<=-2*pi+anglNodes(k)+pi/6)));
      else
         anglDstr(k) = length(find(inlierAngl>=anglNodes(k)-pi/6 & ...
            inlierAngl<=anglNodes(k)+pi/6));
      end
   end

   %Identify the angles that have the minimum population. If the left most of
   % such an angle is to the right of the angle that has the maximum
   % population, we choose the median of such angles to divide the angle
   % curcle so that clusters are always inside the range.
   minVal = min(anglDstr);
   [maxVal,maxInd] = max(anglDstr);
   minInd = sort(find(anglDstr==minVal));
   if minInd(1) > maxInd
      %Choose the closest index from the left to the median.
      dividInd = median(minInd);
      dividInd = max(minInd(find(minInd<=dividInd)));
      dividAngl = anglNodes(dividInd);
      tmp = anglDstr;
      anglDstr(1:end-dividInd+1)   = tmp(dividInd:end);
      anglDstr(end-dividInd+2:end) = tmp(1:dividInd-1);
      tmp = anglNodes;
      anglNodes(1:end-dividInd+1)   = tmp(dividInd:end);
      anglNodes(end-dividInd+2:end) = tmp(1:dividInd-1)+2*pi;

      %Shift 'angl' accordingly.
      ind = find(b.angl<dividAngl);
      b.angl(ind) = b.angl(ind)+2*pi;

      ind = find(inlierAngl<dividAngl);
      inlierAngl(ind) = inlierAngl(ind)+2*pi;
   end

   %Find all local maximum.
   locMaxId = [];
   for k = 2:length(anglNodes)-1
      if (anglDstr(k) > anglDstr(k-1) & anglDstr(k) >= anglDstr(k+1)) | ...
         (anglDstr(k) >= anglDstr(k-1) & anglDstr(k) > anglDstr(k+1))
         if ~isempty(locMaxId) & ...
            min(anglDstr(locMaxId(end)+1:k))>=anglDstr(locMaxId(end))
            if anglDstr(k) > anglDstr(locMaxId(end))
               locMaxId(end) = k;
            else
               locMaxId(end) = floor((locMaxId(end)+k)/2);
            end
         else
            locMaxId = [locMaxId k];
         end
      end
   end

   %Assign angle bin so that its size are big enough (pi/3) and and its center 
   % is local maximum.
   j  = 1;
   jj = 1;
   anglBin(1)  = anglNodes(1);
   selLocMaxId = [];
   anglBinCtr  = [];
   for k = 1:length(locMaxId)
      %To make sure the angle bin size is ~pi/3, we check the distance from
      % current local maximum and the previous angle bin edge.
      if anglNodes(locMaxId(k)) > anglBin(j)+pi/6
         selLocMaxId = [selLocMaxId locMaxId(k)];
         anglBinCtr  = [anglBinCtr anglNodes(locMaxId(k))];

         anglBin(j+1) = anglBinCtr(end)-pi/6;
         anglBin(j+2) = anglBinCtr(end)+pi/6;
         j  = j+2;
      elseif anglNodes(locMaxId(k)) > anglBin(j) & jj > 1
         selLocMaxId = [selLocMaxId locMaxId(k)];
         anglBinCtr  = [anglBinCtr anglNodes(locMaxId(k))];

         anglBin(j)   = (angBinCtr(end)+anglBinCtr(end-1))/2;
         anglBin(j+1) = anglNodes(locMaxId(k))+pi/6;
         j  = j+1;
      elseif ~isempty(selLocMaxId) & ...
         anglDstr(locMaxId(k)) > anglDstr(selLocMaxId(end))
         selLocMaxId(end) = locMaxId(k);
         anglBinCtr(end)  = anglNodes(locMaxId(k));
         anglBin(j-1) = anglBinCtr(end)-pi/6;
         anglBin(j)   = anglBinCtr(end)+pi/6;
      end
   end
   anglBin(end+1) = anglNodes(end);

   j = 1;
   while j < length(anglBin)
      if abs(anglBin(j+1)-anglBin(j)) <= pi/3
         if j == 1
            anglBin(j+1) = [];
         else
            %If two neighbor nodes are close, merge them to the middle point.
            anglBin(j) = (anglBin(j)+anglBin(j+1))/2;
            anglBin(j+1) = [];
         end
         j = j+1;
      else
         %Add nodes.
         nAddedNodes = floor(abs(anglBin(j+1)-anglBin(j))/(pi/3))-1;
         if nAddedNodes > 0
            addedNodes = linspace(anglBin(j),anglBin(j+1),nAddedNodes+2);
            anglBin = [anglBin(1:j) addedNodes(2:end-1) anglBin(j+1:end)];
         end
         j = j+nAddedNodes+2;
      end
   end

   %anglBin  = linspace(-pi,pi,13).';
   nAnglBin = length(anglBin)-1;
   binSet.blockID   = cell(nAnglBin,1);
   binSet.imgPixInd = cell(nAnglBin,1);
   binSet.vecInd    = cell(nAnglBin,1);
   
   vecSet.numObjects = 0;
   img(:)   = 0;
   labelImg = img;
   anglImg  = img;
   for k = 1:nAnglBin
      %Get the index of blocks that belongs to each angle bin.
      %binSet.blockID{k} = b.inlier(find(b.angl(b.inlier)>=anglBin(k) & ...
      %   b.angl(b.inlier)<anglBin(k+1)));
      binSet.blockID{k} = b.inlier(find(inlierAngl>=anglBin(k) & ...
         inlierAngl<anglBin(k+1)));

      binSet.imgPixInd{k} = [];
      if ~isempty(binSet.blockID{k})
         for j = 1:length(binSet.blockID{k})
            blockID = binSet.blockID{k}(j);
            binSet.imgPixInd{k} = [binSet.imgPixInd{k} ...
               b.imgPixInd{blockID}];
         end
      end

      %For each angle bin set population, identify connected clusters
      % (objects) and then put all the objects from all angle bin sets
      % together.
      if ~isempty(binSet.imgPixInd{k})
         %Assign an intensity by the angle
         % bin id for the image.
         img(binSet.imgPixInd{k}) = k;
         anglImg = anglImg + img;
         [L, numObjects] = bwlabel(img);

         nzInd = find(L~=0);
         L(nzInd) = L(nzInd) + vecSet.numObjects;
         labelImg = labelImg + L;

         %Identify the object ID (vecSetID) each block belong to.
         for j = 1:length(binSet.blockID{k})
            blockID = binSet.blockID{k}(j);
            b.vecSetID(blockID) = ...
               labelImg(b.imgPixInd{blockID}(1));
         end
         vecSet.numObjects = vecSet.numObjects+numObjects;

         img(:) = 0;
      end
   end

   if vecSet.numObjects == 0
      return;
   end
   
   %Get the center of each region.
   props = regionprops(labelImg,{'centroid'});
   center = [props.Centroid];
   vecSet.center(:,1) = center(1:2:end);
   vecSet.center(:,2) = center(2:2:end);

   %Identify vector set for each connected cluster (label).
   vecSet.ind          = cell(vecSet.numObjects,1);
   vecSet.avgVc        = NaN*ones(vecSet.numObjects,2);
   vecSet.unitAvgVc    = NaN*ones(vecSet.numObjects,2);
   vecSet.rotUnitAvgVc = NaN*ones(vecSet.numObjects,2);

   vecSet.inlier = [];
   for k = 1:vecSet.numObjects
      vecSet.ind{k} = find(labelImg(indP)==k);

      if ~isempty(vecSet.ind{k})
         %Bundle all the vector indices selected.
         vecSet.inlier = [vecSet.inlier; vecSet.ind{k}];
         
         %'vecSet.avgRotVc' is the vector we use to calcluated the rotation.
         vecSet.avgRotVc(k,:) = mean(rotVc(vecSet.ind{k},:));
         unitAvgVc = vecSet.avgRotVc(k,:)/norm(vecSet.avgRotVc(k,:));

         rotM = [unitAvgVc(2) -unitAvgVc(1);unitAvgVc];

         rotVa(vecSet.ind{k},:) = rotVa(vecSet.ind{k},:)*rotM.';
         rotVh(vecSet.ind{k},:) = rotVh(vecSet.ind{k},:)*rotM.';
         
         vecSet.unitAvgVc(k,:) = unitAvgVc;

         %To check the correctness of the rotation, we also rotate 'vecSet.Vc' 
         % to see if they all align with the vector [0 1].
         vecSet.rotUnitAvgVc(k,:) = vecSet.unitAvgVc(k,:)*rotM.';
      end
   end

   rotVa = rotVa(vecSet.inlier,:);
   rotVh = rotVh(vecSet.inlier,:);
   rotP  = P(vecSet.inlier,:);
   if ~isempty(figH)
      figure(figH(1)); hold off;
      %Blocks that are not used in statistic analysis in shown in black.
      % To stretch the contrast,
      %anglImg(find(anglImg==0)) = -floor(max(anglImg(:))/2);
      %labelImg(find(labelImg==0)) = -floor(max(labelImg(:))/2);

      %stretch the contrast.
      labelImg = labelImg*4;
      labelImg(1,1) = 1.5*max(labelImg(:));

      imshow(labelImg,[]);
      hold on;

      %Plot grid lines
      for k = 1:length(y)
         plot(x-0.5,(y(k)-0.5)*ones(size(x)),'w')
      end
      for k = 1:length(x)
         plot((x(k)-0.5)*ones(size(y)),y-0.5,'w');
      end

      %quiver(P(:,1),P(:,2),rotVa(:,1)*dispScale2,rotVa(:,2)*dispScale2,0,'y');
      %quiver(P(:,1),P(:,2),rotVh(:,1)*dispScale2,rotVh(:,2)*dispScale2,0,'g');
      quiver(P(:,1),P(:,2),Va(:,1)*dispScale,Va(:,2)*dispScale,0,'r');
      quiver(P(:,1),P(:,2),Vh(:,1)*dispScale,Vh(:,2)*dispScale,0,'g');

      quiver(b.center(b.inlier,1),b.center(b.inlier,2), ...
         b.avgVc(b.inlier,1)*dispScale, ...
         b.avgVc(b.inlier,2)*dispScale,0,'w');
      plot(b.center(b.inlier,1),b.center(b.inlier,2),'wo');
      plot(b.center(b.outlier,1),b.center(b.outlier,2),'wo');
      %quiver(b.center(:,1),b.center(:,2),b.avgRotVc(:,1)*corLen, ...
      %   b.avgRotVc(:,2)*corLen,0,'b');

      if length(figH) > 1
         figure(figH(2)); hold off;
         imshow(labelImg,[]); hold on; colormap('jet');

         quiver(rotP(:,1),rotP(:,2),rotVa(:,1)*dispScale,rotVa(:,2)*dispScale,0,'r');
         quiver(rotP(:,1),rotP(:,2),rotVh(:,1)*dispScale,rotVh(:,2)*dispScale,0,'g');
         
         quiver(vecSet.center(:,1),vecSet.center(:,2), ...
            vecSet.unitAvgVc(:,1)*corLen,vecSet.unitAvgVc(:,2)*corLen,0,'w');
         quiver(vecSet.center(:,1),vecSet.center(:,2), ...
            vecSet.rotUnitAvgVc(:,1)*corLen, vecSet.rotUnitAvgVc(:,2)*corLen,0,'k');
      end
   end
   
end

%%% Assemble Data.
avgRotSpdA = mean(sqrt(sum(rotVa.^2,2)));
avgRotSpdH = mean(sqrt(sum(rotVh.^2,2)));
meanVa     = mean(rotVa);
meanVh     = mean(rotVh);
meanVaP    = [meanVa(2) -meanVa(1)];

Da = norm(meanVa);
Dh = norm(meanVh);

alpha = Dh/Da;
dAngl = sign(dot(meanVh,meanVaP))*acos((dot(meanVh,meanVa))/Dh/Da)*180/pi;
S.alpha = alpha;
S.dAngl = dAngl;

