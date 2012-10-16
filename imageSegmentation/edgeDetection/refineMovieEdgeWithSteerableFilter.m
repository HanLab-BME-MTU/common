function prctileUsed = refineMovieEdgeWithSteerableFilter(MD,threshParam,gapCloseParam,doPlot)
%REFINEMOVIEEDGEWITHSTEERABLEFILTER replaces simple masks with masks refined using steerable line filtering
%
%SYNOPSIS refineMovieEdgeWithSteerableFilter(MD,threshParam,gapCloseParam,doPlot)
%
%INPUT  MD    : The movieData object as output by the cell masking and
%               windowing software. Before calling this code,
%               thresholdMovie and refineMovieMask must have been
%               processed. Otherwise the code will crash.
%       threshParam  : Structure with parameters for gradient thresholding:
%           .filterSigma    : Standard deviation for filtering.
%                             Optional. Default: 1.5.
%           .gradPrctile    : Gradient percentile for thresholding.
%                             Optional. Default: [95 90 85 80].
%       gapCloseParam: Structure with parameters for edge gap closing:
%           .maxEdgePairDist: Maximum distance between edge segment pair.
%                             Optional. Default: 5 pixels.
%           .maxThetaDiff   : Maximum angle between gradients of edge
%                             segment pair.
%                             Optional. Default: pi (i.e. full range).
%           .maxC2cAngleThetaDiff: Maximum angle between edge gradient and
%                             perpendicular to centroid-centroid vector.
%                             Optional. Default: pi/2 (i.e. full range).
%           .factorContr    : Contribution of each factor to the edge gap
%                             closing cost. 6 entries for the factors:
%                             (1) distance,
%                             (2) angle between gradients,
%                             (3) angle between gradient and perpendicular
%                                 to centroid-centroid distance,
%                             (4) "edginess" score,
%                             (5) intensity,
%                             (6) lack of asymmetry cost.
%                             Optional. Default: [1 1 1 1 1 1].
%           .edgeType       : Flag indicating edge type:
%                             0 = open edge, i.e. image is of part of a
%                             cell and edge touches image boundary.
%                             1 = closed edge, i.e. image is of whole cell
%                             and segmentation requires finding a closed
%                             contour.
%                             2 = closed edge(s), but of potentially more
%                             than one cell. TO BE IMPLEMENTED IN THE
%                             FUTURE IF NEEDED.
%                             Optional. Default: 0.
%           .fracImageCell  : Fraction of image covered by cell. This
%                             number does not have to be accurate, just
%                             some minimum value to help asses whether
%                             segmentation has been achieved.
%                             Optional. Default: 0.25.
%       doPlot: 1 to plot masks in the end, 2 to also show edge progress,
%               0 to plot nothing. In final plot, refined masks shown in
%               green, original masks shown in blue. Note that this
%               refinement comes on top of the "refineMovieMask"
%               refinement.
%               Optional. Default: 0.
%
%
%OUTPUT prctileUsed: Percentile used for gradient thresholding, one per
%                    frame.
%
%REMARKS The code will copy the original masks and refined masks into new
%        directories called masks_OLD and refined_masks_OLD and will
%        replace the old masks with the ones it produces. After this, one
%        should run refineMovieMasks one more time to get back on track and
%        use the rest of the windowing functions.
%
%Khuloud Jaqaman, November 2011

%% Input/Output

if nargin < 1
    error('refineMovieEdgeWithSteerableFilter: Wrong number of input arguments');
end

%get thresholding parameters, including steerable filter parameters
if nargin < 2 || isempty(threshParam)
    threshParam.filterSigma = 1.5;
    threshParam.gradPrctile = [95 90 85 80];
else
    if ~isfield(threshParam,'fiterSigma')
        threshParam.filterSigma = 1.5;
    end
    if ~isfield(threshParam,'gradPrctile')
        threshParam.gradPrctile = [95 90 85 80];
    end
end

%get edge gap closing parameters
if nargin < 3 || isempty(gapCloseParam)
    gapCloseParam.maxEdgePairDist = 5;
    gapCloseParam.maxThetaDiff = pi;
    gapCloseParam.maxC2cAngleThetaDiff = pi/2;
    gapCloseParam.factorContr = ones(1,6);
    gapCloseParam.edgeType = 0;
    gapCloseParam.fracImageCell = 0.25;
else
    if ~isfield(gapCloseParam,'maxEdgePairDist')
        gapCloseParam.maxEdgePairDist = 5;
    end
    if ~isfield(gapCloseParam,'maxThetaDiff')
        gapCloseParam.maxThetaDiff = pi;
    end
    if ~isfield(gapCloseParam,'maxC2cAngleThetaDiff')
        gapCloseParam.maxC2cAngleThetaDiff = pi/2;
    end
    if ~isfield(gapCloseParam,'factorContr')
        gapCloseParam.factorContr = ones(1,6);
    end
    if ~isfield(gapCloseParam,'edgeType')
        gapCloseParam.edgeType = 0;
    end
    if ~isfield(gapCloseParam,'fracImageCell')
        gapCloseParam.fracImageCell = 0.25;
    end
end

%check whether/what to plot
if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end

%get image and analysis directories
imageDir = MD.channels_.channelPath_;
analysisDir = MD.movieDataPath_;

%make new directories and copy old masks to them
masksDir = [analysisDir filesep 'masks'];
masksOldDir = [analysisDir filesep 'masks_OLD'];
mkdir(masksOldDir);
copyfile([masksDir filesep '*'],masksOldDir);
refinedMasksDir = [analysisDir filesep 'refined_masks'];
refinedMasksOldDir = [analysisDir filesep 'refined_masks_OLD'];
mkdir(refinedMasksOldDir);
copyfile([refinedMasksDir filesep '*'],refinedMasksOldDir);

%get image and mask file listings
imageFileListing = dir([imageDir filesep '*.tif']);
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.tiff']);
end
masksDirFull = [masksDir filesep 'masks_for_channel_1'];
maskFileListing = dir([masksDirFull filesep '*.tif']);
if isempty(maskFileListing)
    maskFileListing = dir([masksDirFull filesep '*.tiff']);
end
refinedMasksDirFull = [refinedMasksDir filesep 'refined_masks_for_channel_1'];
refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tif']);
if isempty(refinedMaskFileListing)
    refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tiff']);
end

%get number of files
numFiles = length(imageFileListing);

%% Mask refinement

wtBar = waitbar(0,'Please wait, refining edge with steerable filter ...'); 

%refine masks using steerable line filter
prctileUsed = NaN(numFiles,1);
for iFile = 1 : numFiles
    
    waitbar(iFile/numFiles,wtBar,'Please wait, refining edge with steerable filter ...'); 
    
    %read image and mask
    image = double(imread(fullfile(imageDir,imageFileListing(iFile).name)));
    mask0 = double(imread(fullfile(refinedMasksDirFull,refinedMaskFileListing(iFile).name)));
    
    %call refinement function
    [mask,prctileUsed(iFile)] = refineEdgeWithSteerableFilter(mask0,image,...
        threshParam,gapCloseParam,doPlot);
    
    %store new mask
    imwrite(mask,fullfile(masksDirFull,maskFileListing(iFile).name),'tif');
    
    if prctileUsed(iFile)==-1
        disp(['bad mask Image ' num2str(iFile)]);
    end
    
end

close(wtBar)

%% ~~~ the end ~~~

