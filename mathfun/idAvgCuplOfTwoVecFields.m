function S = idAvgCuplOfTwoVecFields(Va,Vh,P,varargin)
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
%    'rawVF'          : A cell array of the original two vector fields used 
%                       to get the interpolated 'Va' and 'Vh'. Vector field
%                       is given in the format [x y vx vy]. 
%                       Note: It makes more sense to use the original vector
%                       population for speed, coherence and population test.
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
%        'aCohrS'    : The cohrence score of 'Va'.
%        'hCohrS'    : The cohrence score of 'Vh'.
%        'cohrPPLP'  : The percentage of population that pass the coherence
%                      test.
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
%    below 'smSpdCohrThreshold'.

%Default paramters.
bSize               = 30;
scaling             = 'on';
smSpdThreshold      = 0;
smSpdCohrThreshold  = 0.05;
bigSpdCohrThreshold = 0.05;
pplThreshold        = 3;
figH                = [];
img                 = [];
Vc                  = Va;
dispScale           = 15;

%Default return.
S.alpha    = NaN;
S.dAngl    = NaN;
S.aCohrS   = NaN;
S.hCohrS   = NaN;
S.cohrPPLP = 0;

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

%Default raw vector fields.
rawPa = P;
rawVa = Va;
rawPh = P;
rawVh = Vh;

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
      case 'rawVF'
         rawPa = varargin{k+1}{1}(:,1:2);
         rawVa = varargin{k+1}{1}(:,3:4);
         rawPh = varargin{k+1}{2}(:,1:2);
         rawVh = varargin{k+1}{2}(:,3:4);

         %Remove NaN from rawVa and rawVh.
         nanInd = find(isnan(rawVa(:,1)) | isnan(rawVa(:,2)));
         rawPa(nanInd,:) = [];
         rawVa(nanInd,:) = [];

         nanInd = find(isnan(rawVh(:,1)) | isnan(rawVh(:,2)));
         rawPh(nanInd,:) = [];
         rawVh(nanInd,:) = [];
      case 'smSpdThreshold'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            smSpdThreshold = varargin{k+1};
         end
      case 'cohrThreshold'
         if ~isempty(varargin{k+1}) & ~isnan(varargin{k+1})
            if length(varargin{k+1}) == 2
               smSpdCohrThreshold = varargin{k+1}(1);
               bigSpdCohrThreshold = varargin{k+1}(2);
            else
               smSpdCohrThreshold = varargin{k+1};
               bigSpdCohrThreshold = varargin{k+1};
            end
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
         %Remove NaN.
         nanInd = find(isnan(Vc(:,1)) | isnan(Vc(:,2)));
         Vc(nanInd,:) = [];
         if isempty(Vc)
            return;
         end
      otherwise
         error(['Optional parameter '' varargin{k} '' is not recogonized.']);
   end
end

%Get the raw population that are in the neighborhood of P.
minX = min(P(:,1));
maxX = max(P(:,1));
minY = min(P(:,2));
maxY = max(P(:,2));

if isinf(bSize) | isnan(bSize)
   infLen = 5;
else
   infLen = bSize/2;
end
rawAInd = find(rawPa(:,1)<maxX+infLen & ...
   rawPa(:,1)>minX-infLen & ...
   rawPa(:,2)<maxY+infLen & ...
   rawPa(:,2)>minY-infLen); 

if isempty(rawAInd)
   return;
end
rawPa = rawPa(rawAInd,:);
rawVa = rawVa(rawAInd,:);

rawHInd = find(rawPh(:,1)<maxX+infLen & ...
   rawPh(:,1)>minX-infLen & ...
   rawPh(:,2)<maxY+infLen & ...
   rawPh(:,2)>minY-infLen); 

if isempty(rawHInd)
   return;
end
rawPh = rawPh(rawHInd,:);
rawVh = rawVh(rawHInd,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation to preserve spatial heterogeniety.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create an image area that covers the vector field.
corLen = ceil(bSize/2);

cSpd = sqrt(sum(Vc.^2,2));
aSpd = sqrt(sum(Va.^2,2));
hSpd = sqrt(sum(Vh.^2,2));

rawASpd = sqrt(sum(rawVa.^2,2));
rawHSpd = sqrt(sum(rawVh.^2,2));

%To avoid dividing by zero.
smNumber = min([1e-6 1e-6*[median(aSpd) median(hSpd) ...
   median(rawASpd) median(rawHSpd)]]);
if smNumber == 0
   smNumber = 1e-6;
end

unitRawVa = NaN*ones(size(rawVa));
unitRawVh = NaN*ones(size(rawVh));

nzInd = find(rawASpd>=smNumber);
if ~isempty(nzInd)
   unitRawVa(nzInd,1) = rawVa(nzInd,1)./rawASpd(nzInd);
   unitRawVa(nzInd,2) = rawVa(nzInd,2)./rawASpd(nzInd);
end

nzInd = find(rawHSpd>=smNumber);
if ~isempty(nzInd)
   unitRawVh(nzInd,1) = rawVh(nzInd,1)./rawHSpd(nzInd);
   unitRawVh(nzInd,2) = rawVh(nzInd,2)./rawHSpd(nzInd);
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

if isinf(corLen) | (maxX-minX <= bSize & maxY-minY <= bSize) 
   %When corLen is Inf or the whole population is within the size of a
   % block, we calculate the average of the whole population. 
   % We essentially do nothing here except for the speed, coherence
   % and population test.

   avgRawASpd = mean(rawASpd);

   numRawVa = find(~isnan(unitRawVa(:,1)) & ~isnan(unitRawVa(:,2)));
   if ~isempty(numRawVa)
      S.aCohrS = norm(mean(unitRawVa(numRawVa,:),1));
   end

   numRawVh = find(~isnan(unitRawVh(:,1)) & ~isnan(unitRawVh(:,2)));
   if ~isempty(numRawVh)
      S.hCohrS = norm(mean(unitRawVh(numRawVh,:),1));
   end

   %Speed, coherence and population test.
   if isnan(S.aCohrS) | ...
      length(numRawVa) <= pplThreshold | length(numRawVh) <= pplThreshold | ...
      (avgRawASpd > smSpdThreshold & S.aCohrS <= bigSpdCohrThreshold) | ...
      (avgRawASpd <= smSpdThreshold & S.aCohrS <= smSpdCohrThreshold)
      return;
   end

   S.cohrPPLP = 1;
else
   if isempty(img)
      x = minX-corLen:2*corLen:maxX+3*corLen;
      y = minY-corLen:2*corLen:maxY+3*corLen;

      %Translate so that the image area starts with pixel (1,1).
      P(:,1) = P(:,1)-x(1)+1;
      P(:,2) = P(:,2)-y(1)+1;

      rawPa(:,1) = rawPa(:,1)-x(1)+1;
      rawPa(:,2) = rawPa(:,2)-y(1)+1;
      rawPh(:,1) = rawPh(:,1)-x(1)+1;
      rawPh(:,2) = rawPh(:,2)-y(1)+1;

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
   [imgHeight,imgWidth] = size(img);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Divide the image area into blocks and identify the vector population for
   % each block.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   nGridY  = length(y);
   nGridX  = length(x);
   nBlocks = (nGridX-1)*(nGridY-1);

   %To cover all possible clusters of points, we have 4 layers of blocks. Each
   % of them is a half block shift of the base layer.
   %Block shift for each layer. Format: [y x].
   bShift  = bSize*[0 0; ...
                    0 0.5; ...
                    0.5 0; ...
                    0.5 0.5];

%    bShift = [0 0];
   nLayers = size(bShift,1);

   %Round P, rawPa and rawPh to integer pixel and get their linear index
   % relative to 'img'.
   rndP = round(P);
   indP = sub2ind(size(img),rndP(:,2),rndP(:,1));

   rndPa = round(rawPa);
   rndPh = round(rawPh);
   indPa = sub2ind(size(img),rndPa(:,2),rndPa(:,1));
   indPh = sub2ind(size(img),rndPh(:,2),rndPh(:,1));

   %We need the image index of each block in each layer later.
   b.imgPixInd = cell(nBlocks,nLayers);
   k = 0;
   for jx = 1:nGridX-1
      for jy = 1:nGridY-1
         k = k+1;
         for jl = 1:nLayers
            [bX bY] = meshgrid([x(jx):x(jx+1)-1]+bShift(jl,2), ...
               [y(jy):y(jy+1)-1]+bShift(jl,1));

            %Make sure we only get the part of block that is inside the iamge.
            in = find(bX>=1 & bX<=imgWidth & bY>=1 & bY<=imgHeight);
            b.imgPixInd{k,jl} = sub2ind(size(img),bY(in),bX(in));
         end
      end
   end

   %Index of vectors inside each block of each layer.
   b.vecInd     = cell(nBlocks,nLayers); %Vc
   b.rVaInd     = cell(nBlocks,nLayers); %rawVa
   b.rVhInd     = cell(nBlocks,nLayers); %rawVh

   %Center of each block in the base layer.
   b.center     = zeros(nBlocks,2);

   %Average flow angle, coupling vector, speed of raw vector field.
   b.angl       = NaN*ones(nBlocks,nLayers);
   b.avgVc      = NaN*ones(nBlocks,2,nLayers);
   b.avgRawASpd = NaN*ones(nBlocks,nLayers);

   %Coherence of the two flow channels (Va and Vh) in each block of each layer.
   b.aCohrS     = NaN*ones(nBlocks,nLayers);
   b.hCohrS     = NaN*ones(nBlocks,nLayers);

   %The index of outlier blocks (do not pass test) and inlier blocks (pass
   % test) at each layer.
   b.outlier    = cell(1,nLayers);
   b.inlier     = cell(1,nLayers);
   for jl = 1:nLayers
      b.outlier{jl} = [];
      b.inlier{jl}  = [];
   end

   %Identify vector index, calculate coherence score and do coherence test for
   % each block population.
   for k = 1:nBlocks
      %Identify the starting index and the center of the block in the 
      % base layer.
      [yI xI]     = ind2sub([nGridY-1,nGridX-1],k);
      b.center(k,:) = [(x(xI)+x(xI+1))/2 (y(yI)+y(yI+1))/2];

      for jl = 1:nLayers
         %Identify the index of vectors that are inside each block.
         %First, block boundary in y and x direction.
         yL = y(yI)+bShift(jl,1);
         yR = y(yI+1)+bShift(jl,1);
         xL = x(xI)+bShift(jl,2);
         xR = x(xI+1)+bShift(jl,2);
         b.vecInd{k,jl} = find(rndP(:,2)>=yL & rndP(:,2) < yR & ...
            rndP(:,1)>=xL & rndP(:,1)<xR);

         b.rVaInd{k,jl} = find(rndPa(:,2)>=yL & rndPa(:,2) < yR & ...
            rndPa(:,1)>=xL & rndPa(:,1)<xR);
         b.rVhInd{k,jl} = find(rndPh(:,2)>=yL & rndPh(:,2) < yR & ...
            rndPh(:,1)>=xL & rndPh(:,1)<xR);

         if ~isempty(b.vecInd{k,jl})
            %Calculate the angle of the average vector for the block.
            % To avoid dividing by zero, we use 'smNumber' defined above.
            b.avgVc(k,:,jl) = mean(Vc(b.vecInd{k,jl},:),1);
            b.angl(k,jl)    = sign(b.avgVc(k,2,jl))*acos(b.avgVc(k,1,jl)/ ...
                              max(smNumber,norm(b.avgVc(k,:,jl))));

            %Calculate the coherence score for each block.
            if ~isempty(b.rVaInd{k,jl})
               nInd = b.rVaInd{k,jl}(find(~isnan(unitRawVa(b.rVaInd{k,jl},1)) & ...
                  ~isnan(unitRawVa(b.rVaInd{k,jl},2))));
               if ~isempty(nInd)
                  b.aCohrS(k,jl) = norm(mean(unitRawVa(nInd,:),1));

                  %Coherence score only makes sense when there is enough
                  % population. So, we scale it by the ratio to
                  % 'pplThreshold'. If there is only one vector, the score is
                  % set to zero. If there is more than 'pplThreshold' vectors,
                  % it is not scaled (multiplied by 1).
                  b.aCohrS(k,jl) = b.aCohrS(k,jl)/pplThreshold* ...
                     min(pplThreshold,length(nInd)-1);
               end
               b.avgRawASpd(k,jl) = mean(rawASpd(b.rVaInd{k,jl}));
            end

            if ~isempty(b.rVhInd{k,jl})
               nInd = b.rVhInd{k,jl}(find(~isnan(unitRawVh(b.rVhInd{k,jl},1)) & ...
                  ~isnan(unitRawVh(b.rVhInd{k,jl},2))));
               if ~isempty(nInd)
                  b.hCohrS(k,jl) = norm(mean(unitRawVh(nInd,:),1));
                  b.hCohrS(k,jl) = b.hCohrS(k,jl)/pplThreshold* ...
                     min(pplThreshold,length(nInd)-1);
               end
            end

            %Speed, coherence and population test.
            if isnan(b.aCohrS(k,jl)) | ...
               length(b.rVaInd{k,jl}) <= pplThreshold | ...
               length(b.rVhInd{k,jl}) <= pplThreshold | ...
               (b.avgRawASpd(k,jl) > smSpdThreshold & b.aCohrS(k,jl) <= bigSpdCohrThreshold) | ...
               (b.avgRawASpd(k,jl) <= smSpdThreshold & b.aCohrS(k,jl) <= smSpdCohrThreshold)
               b.outlier{jl} = [b.outlier{jl}; k];
            else
               b.inlier{jl} = [b.inlier{jl}; k];
            end
         end
      end
   end
   
   %Calculate the mean coherence score.
   numberInd = find(~isnan(b.aCohrS));
   if ~isempty(numberInd)
      S.aCohrS = mean(b.aCohrS(numberInd));
   end
   numberInd = find(~isnan(b.hCohrS));
   if ~isempty(numberInd)
      S.hCohrS = mean(b.hCohrS(numberInd));
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %A smart angle bin assignment:
   % Group blocks into a set of angle bins. We assign angle bins according the
   % natural angle distribution. This is a little bit tricky to code.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %Identify clusters of angle distribution.
   % First, get the population that are withing pi/6 (30 degree) from each
   % angle node in 'anglNodes' which is a vector of evenly distributed angle 
   % points between -pi and pi. 
   anglNodes = linspace(-pi,pi,73);
   anglDstr  = zeros(size(anglNodes));
   
   %Combine the angle of inliner blocks from all layers.
   inlierAngl = [];
   for jl = 1:nLayers
      inlierAngl = [inlierAngl; b.angl(b.inlier{jl},jl)];
   end
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

   %Merge close nodes or add nodes to too big angle bin.
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
   %%% end of angle bin assignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   nAnglBin = length(anglBin)-1;

   %Index of blocks whose angle belongs to each angle bin in each layer.
   binSet.blockID   = cell(nAnglBin,nLayers);

   %The image index under the blocks that belong to each bin.
   binSet.imgPixInd = cell(nAnglBin,1);

   %Index of vectors that belong to each bin.
   binSet.vecInd    = cell(nAnglBin,nLayers);
   
   %Make sure intially 'img = 0'.
   img(:)        = 0;

   %To solve the problem of overlaying blocks, we assign the coherence score
   % of each block to the pixels it covers and compare the score between
   % layers.
   %lastACohrSImg = img;
   %thisACohrSImg = img;
   aCohrSImg = img;
   for k = 1:nAnglBin
      binSet.imgPixInd{k} = [];
      for jl = 1:nLayers
         %Get the index of blocks that belongs to each angle bin.
         binSet.blockID{k,jl} = b.inlier{jl}( ...
            find(b.angl(b.inlier{jl},jl)>=anglBin(k) & ...
            b.angl(b.inlier{jl},jl)<anglBin(k+1)));

         %Get the image pixels covered by all the blocks (from all layers) that
         % belongs to each angle bin.
         if ~isempty(binSet.blockID{k,jl})
            for j = 1:length(binSet.blockID{k,jl})
               blockID = binSet.blockID{k,jl}(j);
               bImgPix = b.imgPixInd{blockID,jl};
               binSet.imgPixInd{k} = [binSet.imgPixInd{k}; bImgPix];

               %We only update the cohrence score for pixels whose current
               % score is less than the coherence score of the block at this
               % layer that is going to cover it.
               updInd = bImgPix(find(aCohrSImg(bImgPix) < ...
                  b.aCohrS(blockID,jl)));
               %thisACohrSImg(updInd) = b.aCohrS(blockID,jl);
               aCohrSImg(updInd) = b.aCohrS(blockID,jl);
               img(updInd) = k;
            end
         end
      end

%      bsImgPix = binSet.imgPixInd{k};
%      if ~isempty(bsImgPix)
%         %Assign an intensity by the angle bin id for these image pixles.
%         updInd = bsImgPix( ...
%            find(thisACohrSImg(bsImgPix)>lastACohrSImg(bsImgPix)));
%         img(updInd) = k;
%      end
%      lastACohrSImg = thisACohrSImg;
   end

   %Identify clusters of connected blocks that belongs to the same angle bin.
   % Each of such clusters define a subpopulation of vectors.
   nVecSets = 0;
   labelImg = zeros(size(img));
   for k = 1:nAnglBin
      %For each angle bin set population, identify connected block 
      % clusters (objects).
      L = img;
      L(find(img~=k)) = 0;
      [L, numObjects] = bwlabel(L);

      nzInd = find(L~=0);
      L(nzInd) = L(nzInd) + nVecSets;

      %Then, put all the objects from all angle bin sets
      % together. This is represented by 'labelImg' which is all 
      % we need to identify clusters of vectors whose angles 
      % are close.
      labelImg = labelImg + L;

      nVecSets = nVecSets+numObjects;
   end

   if nVecSets == 0
      return;
   end
   
   %Get the center of each region (cluster of vectors).
   props = regionprops(labelImg,{'centroid'});
   center = [props.Centroid];
   vecSet.center(:,1) = center(1:2:end);
   vecSet.center(:,2) = center(2:2:end);

   %The index of vectors that belong to each connected cluster (label).
   vecSet.ind          = cell(nVecSets,1);
   vecSet.avgRotVc     = NaN*ones(nVecSets,2);
   vecSet.unitAvgVc    = NaN*ones(nVecSets,2);
   vecSet.rotUnitAvgVc = NaN*ones(nVecSets,2);

   vecSet.inlier = [];
   for k = 1:nVecSets
      vecSet.ind{k} = find(labelImg(indP)==k);

      if length(vecSet.ind{k}) > pplThreshold
         %We only consider vector sets whose population is above the
         % population threshold, 'pplThreshold'.
         %Bundle all the vector indices selected.
         vecSet.inlier = [vecSet.inlier; vecSet.ind{k}];
         
         %'vecSet.avgRotVc' is the vector we use to calcluated the rotation.
         vecSet.avgRotVc(k,:) = mean(rotVc(vecSet.ind{k},:),1);
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

   if ~isempty(figH)
      vecSet.rawVaInlier = find(labelImg(indPa)~=0);
      vecSet.rawVhInlier = find(labelImg(indPh)~=0);

      figure(figH(1)); hold off;
      %Blocks that are not used in statistic analysis in shown in black.
      % We stretch the contrast so that each angle bin cluster can be
      % distinguished.
      minV    = min(img(:));
      maxV    = max(img(:));
      imgDisp = img;

      nzInd = find(img~=0);
      imgDisp(nzInd) = (img(nzInd)-2*minV+maxV);
      imgDisp(1,1)   = 3*(maxV-minV);

      imshow(4*imgDisp,[]);
      hold on;

      %Plot grid lines
      for k = 1:length(y)
         plot(x-0.5,(y(k)-0.5)*ones(size(x)),'w')
      end
      for k = 1:length(x)
         plot((x(k)-0.5)*ones(size(y)),y-0.5,'w');
      end

      quiver(rawPa(:,1),rawPa(:,2),rawVa(:,1)*dispScale, ...
         rawVa(:,2)*dispScale,0,'y');
      quiver(rawPh(:,1),rawPh(:,2),rawVh(:,1)*dispScale, ...
         rawVh(:,2)*dispScale,0,'b');

      inlier = vecSet.rawVaInlier;
      quiver(rawPa(inlier,1),rawPa(inlier,2),rawVa(inlier,1)*dispScale, ...
         rawVa(inlier,2)*dispScale,0,'r');
      inlier = vecSet.rawVhInlier;
      quiver(rawPh(inlier,1),rawPh(inlier,2),rawVh(inlier,1)*dispScale, ...
         rawVh(inlier,2)*dispScale,0,'g');

      quiver(b.center(b.inlier{1},1),b.center(b.inlier{1},2), ...
         b.avgVc(b.inlier{1},1,1)*dispScale, ...
         b.avgVc(b.inlier{1},2,1)*dispScale,0,'w');
      plot(b.center(b.inlier{1},1),b.center(b.inlier{1},2),'wo');
      plot(b.center(b.outlier{1},1),b.center(b.outlier{1},2),'ro');
      tH = text(b.center(b.outlier{1},1),b.center(b.outlier{1},2)-5,num2str(b.outlier{1}));
      set(tH,'color','r');

      if length(figH) > 1
         %stretch the contrast.
         minV      = min(labelImg(:));
         maxV      = max(labelImg(:));
         labelDisp = labelImg;

         nzInd = find(labelImg~=0);
         labelDisp(nzInd) = labelDisp(nzInd)-2*minV+maxV;
         labelDisp(1,1)   = 3*(maxV-minV);

         figure(figH(2)); hold off;
         imshow(4*labelDisp,[]); hold on; colormap('jet');

         quiver(P(:,1),P(:,2),rotVa(:,1)*dispScale,rotVa(:,2)*dispScale,0,'r');
         quiver(P(:,1),P(:,2),rotVh(:,1)*dispScale,rotVh(:,2)*dispScale,0,'g');
         
         inlier = vecSet.inlier;
         quiver(P(inlier,1),P(inlier,2),rotVa(inlier,1)*dispScale, ...
            rotVa(inlier,2)*dispScale,0,'w');
         quiver(P(inlier,1),P(inlier,2),rotVh(inlier,1)*dispScale, ...
            rotVh(inlier,2)*dispScale,0,'k');

         quiver(vecSet.center(:,1),vecSet.center(:,2), ...
            vecSet.unitAvgVc(:,1)*corLen,vecSet.unitAvgVc(:,2)*corLen,0,'w');
         quiver(vecSet.center(:,1),vecSet.center(:,2), ...
            vecSet.rotUnitAvgVc(:,1)*corLen, vecSet.rotUnitAvgVc(:,2)*corLen,0,'k');
         plot(vecSet.center(:,1),vecSet.center(:,2),'wo');
      end
   end

   S.cohrPPLP = length(vecSet.inlier)/size(rotVa,1);
   rotVa = rotVa(vecSet.inlier,:);
   rotVh = rotVh(vecSet.inlier,:);
end

%%% Assemble Data.
%avgRotSpdA = mean(sqrt(sum(rotVa.^2,2)));
%avgRotSpdH = mean(sqrt(sum(rotVh.^2,2)));
meanVa     = mean(rotVa,1);
meanVh     = mean(rotVh,1);
meanVaP    = [meanVa(2) -meanVa(1)];

Da = norm(meanVa);
Dh = norm(meanVh);

alpha = Dh/Da;
dAngl = sign(dot(meanVh,meanVaP))*acos((dot(meanVh,meanVa))/Dh/Da)*180/pi;
S.alpha = alpha;
S.dAngl = dAngl;

