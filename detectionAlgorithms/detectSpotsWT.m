function [frameInfo, imgDenoised] = detectSpotsWT(img, iclean, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addOptional('iclean', 0);
ip.addParamValue('LegacyMode', false);
ip.addParamValue('BackwardsCompatible', true, @islogical); % to match output of main283AUTO
ip.parse(img, iclean, varargin{:});

iclean = ip.Results.iclean;

% Wavelet decomposition calculated up to scale k
k = 4;

if ip.Results.LegacyMode
    denoisefun = @wstarck2;
else
    denoisefun = @awtDenoising;    
end

% NOTE: this file is identical to main283AUTO_standalone in function.
% The output vectors elements may be in a different order due to implementation differences.

% Reconstructs a filtered image using iterative filtering from significant
% coefficients ( see Book by Starck et al, p. 66)

%Then calculates a binary image and the product of the binary image and the
%reconstructed image

% It then finds and plots ALL the local maxima on each cluster (at least
% one maximum per cluster)

%If a cluster has more than one maximum, it divides them into primary and
%secondary (The secondaries are listed at the end.)

%It calculates the intensity of each maximum and the total intensity and size (i.e. area) of
%each cluster.

%   INPUT
%             img : input image
%
%  OPTIONS
%         iclean : 0 when only minimal cleaning of very small clusters is necessary (default)
%                  1 when some extra cleaning is necessary
%                  2 when some more extra cleaning is necessary
%                 -1 when even minimal cleaning is not necessary
%
%  OUTPUT
%
%  Cluster-specific values:  
%
%     xav, yav   : centroid coordinates for the clusters
%     num        : number of clusters
%
%  Local maxima (each cluster can have several)
%     xmax, ymax : coordinates of local maxima 
%     inn        : intensity of local maxima
%     intot      : intensity of the cluster (sum of all pixels in cluster)
%                : this is repeated for each maxima within the cluster (change in future version?)
%     csize      : number of pixels in each cluster (same as above)
%     lxm        : number of local maxima in cluster (same again)
%     labl       : cluster to which local maximum belong (i.e., cluster label) 
%     nmax       : number of local maxima in the image
%

% Last modified: 08/26/12 by Francois Aguet

[ny,nx,nz] = size(img);
imgDenoised = zeros(ny,nx,nz, 'uint8');

frameInfo(1:nz) = struct('ymax', [], 'xmax', [], 'inn', [], 'yav', [], 'xav', [],...
    'intot', [], 'csize', [], 'lxm', [], 'labl', [], 'num', [], 'nmax', []);

% loop through frames in stack
for ix = 1:nz
    
    frame = double(img(:,:,ix));
    
    mxo=max(frame(:));
    mno=min(frame(:));
    frame = frame*450/mxo; % to renormalize the maximum intensity to a reasonable value.
       
    rw = denoisefun(frame,k,0); % CALCULATES reconstructed image up to order k without kth detail(irecon=0)
    er = frame - rw;  % This is the error image
    
    % Now do iterative filtering
    delta = 1;
    n=0;
    sig1 = std(er(:));
    while delta > 0.002
        
        n = n+1;
        rw = rw + denoisefun(er,k,0);
        
        er = frame - rw;
        sig2 = std(er(:));
        delta = abs(sig2-sig1)/sig2;
        sig1 = sig2;
    end
    if ip.Results.LegacyMode
        wm = wstarck222(rw,k); % calculate the multiscale product of the details
    else
        W = abs(awt(rw, k));
        wm = prod(W(:,:,1:k-1),3);
    end
    rw = rw-min(rw(:));
    
    % size of mask for local average = (2*bs+1)by(2*bs+1)
    [av,sg] = localAvgStd2D(rw, 9);

    
    %create binary image
    rwb = zeros(ny, nx);
    rwb((rw >= av+0.5*sg) & (rw.*wm >= mean(rw(:)))) = 1;
    
    rwb=bwmorph(rwb,'clean'); %to get rid of isolated pixels (noise)
    rwb=bwmorph(rwb,'fill'); %to fill up empty internal pixels
    rwb=bwmorph(rwb,'thicken'); %to make larger clusters because of the harsh cutoff criteria above
    rwb=bwmorph(rwb,'spur'); % to remove single pixels 8-attached to clusters
    rwb=bwmorph(rwb,'spur');
    rwb=bwmorph(rwb,'clean');% to remove any remaining isolated pixels after applying spur twice
    if iclean > 0
        rwb=bwmorph(rwb,'erode'); %extra cleaning of small spots
        if iclean >=2
            rwb=bwmorph(rwb,'spur'); %extra cleaning of small spots % for extraextra cleaning
        end
        rwb=bwmorph(rwb,'clean');%extra cleaning of small spots
        rwb=bwmorph(rwb,'thicken');%extra cleaning of small spots
    end
    rw = rw*(mxo-mno)/max(rw(:)); % normalize maximum intensity to original (maximum intensity-background intensity).
    rw = rw.*rwb; %Outside the clusters the intensity is set to zero exactly
    
    
    % Connected components -> clusters of significant pixels
    CC = bwconncomp(rwb, 8);
    num = CC.NumObjects;
    labelMat = double(labelmatrix(CC));
    
    % calculate local maxima
    lm = locmax2d(rw, [9 9]); % find the location of the maxima of the clusters
    [ymax, xmax] = find(lm>0); % coordinates of the local maxima
    
    labl = labelMat(sub2ind([ny nx],ymax,xmax));
    
    % clusters for which local maxima detection failed
    idx = setdiff(1:num, labl);
    
    % determine local max coordinates for these clusters
    idx = cellfun(@(i) i(find(rw(i)==max(rw(i)),1,'first')), CC.PixelIdxList(idx));
    [ymax2, xmax2] = ind2sub([ny nx], idx);
    ymax = [ymax; ymax2']; %#ok<AGROW>
    xmax = [xmax; xmax2']; %#ok<AGROW>
    nmax = length(xmax);

    % update labels for all local maxima
    labl = labelMat(sub2ind([ny nx], ymax, xmax));
    
    % retrieve significant pixel clusters and associated properties
    stats = regionprops(CC, rw, 'WeightedCentroid');
    centroids = vertcat(stats.WeightedCentroid);
    xav = centroids(:,1);
    yav = centroids(:,2);
    
    % per-cluster measurements
    % total cluster intensity
    intot = cellfun(@(i) sum(rw(i)), CC.PixelIdxList)';
    % # pixels in cluster
    csize = cellfun(@(i) numel(i), CC.PixelIdxList)';
    
    % # local max. in each cluster
    lxm = getMultiplicity(labl)';
    
    % intensity of each local maximum
    inn = rw(sub2ind([ny nx], ymax, xmax));
    
    % if backwards-compatible output is desired
    if ip.Results.BackwardsCompatible
        intot = intot(labl);
        csize = csize(labl);
        lxm = lxm(labl);
    end
    
    frameInfo(ix).ymax = ymax;
    frameInfo(ix).xmax = xmax;
    frameInfo(ix).inn = inn;
    frameInfo(ix).yav = yav;
    frameInfo(ix).xav = xav;
    frameInfo(ix).intot = intot;
    frameInfo(ix).csize = csize;
    frameInfo(ix).lxm = lxm;
    frameInfo(ix).labl = labl;
    frameInfo(ix).num = num;
    frameInfo(ix).nmax = nmax;
   
    % eliminate the frequent error messages by rescaling to uint8
    rw8 = uint8(255*(rw/max(rw(:))));
    imgDenoised(:,:,ix) = rw8;
end
