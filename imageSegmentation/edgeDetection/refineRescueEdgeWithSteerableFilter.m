function prctileUsed = refineRescueEdgeWithSteerableFilter(MD,frames2rescue,...
    threshParam,gapCloseParam,doPlot,movieInfo,numSPTFrames,meanBkg)
%refineRescueEdgeWithSteerableFilter like refineMovieEdgeWithSteerableFilter but for specific frames only, can be run only after refineMovieEdgeWithSteerableFilte
%
%SYNOPSIS prctileUsed = refineRescueEdgeWithSteerableFilter(MD,frames2rescue,...
%    threshParam,gapCloseParam,doPlot,movieInfo,numSPTFrames,meanBkg)
%
%INPUT  MD    : The movieData object as output by the cell masking and
%               windowing software. Before calling this code,
%               thresholdMovie and refineMovieMask must have been
%               processed. Otherwise the code will crash.
%       frames2rescue: Row vector of indices of frames to rescue. 
%       threshParam  : Structure with parameters for gradient thresholding:
%           .filterSigma    : Standard deviation for filtering.
%                             Optional. Default: 1.5.
%           .gradPrctile    : Gradient percentile for thresholding.
%                             Optional. Default: [95 90 85 80].
%       gapCloseParam: Structure with parameters for edge gap closing:
%           .maxEdgePairDist: Maximum distance between edge segment pair.
%                             Optional. Default: 5 pixels.
%           .factorContr    : Contribution of each factor to the edge gap
%                             closing cost. 6 entries for the factors:
%                             (1) distance,
%                             (2) angle between gradients,
%                             (3) angle between gradient and perpendicular
%                                 to centroid-centroid distance,
%                             (4) "edginess" score,
%                             (5) lack of asymmetry cost.
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
%       movieInfo: Detection output with particle positions.
%                  Optional. If not input, information not used.
%       numSPTFrames: Number of SPT frames between edge frames.
%                  Only needed if movieInfo is also input.
%                  Optional. Default: 400.
%       meanBkg  : Mean background intensity close to the cell edge.
%                  Optional. If not input, information not used.
%
%OUTPUT prctileUsed: Percentile used for gradient thresholding, one per
%                    frame.
%
%REMARKS ...
%
%Khuloud Jaqaman, November 2011

%% Input/Output

if nargin < 2
    error('refineRescueEdgeWithSteerableFilter: Wrong number of input arguments');
end

%get thresholding parameters, including steerable filter parameters
if nargin < 3 || isempty(threshParam)
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
if nargin < 4 || isempty(gapCloseParam)
    gapCloseParam.maxEdgePairDist = 5;
    gapCloseParam.factorContr = ones(1,6);
    gapCloseParam.edgeType = 0;
    gapCloseParam.fracImageCell = 0.25;
else
    if ~isfield(gapCloseParam,'maxEdgePairDist')
        gapCloseParam.maxEdgePairDist = 5;
    end
    if ~isfield(gapCloseParam,'factorContr')
        gapCloseParam.factorContr = ones(1,5);
    end
    if ~isfield(gapCloseParam,'edgeType')
        gapCloseParam.edgeType = 0;
    end
    if ~isfield(gapCloseParam,'fracImageCell')
        gapCloseParam.fracImageCell = 0.25;
    end
end

%check whether/what to plot
if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end

%check particle information
if nargin < 6 || isempty(movieInfo)
    movieInfo = [];
end
if nargin < 7 || isempty(numSPTFrames)
    numSPTFrames = 400;
end

if nargin < 8 || isempty(meanBkg)
    meanBkg = [];
end

%get image and analysis directories
imageDir = MD.channels_.channelPath_;
analysisDir = MD.movieDataPath_;

%make new directories and copy old masks to them
masksDir = [analysisDir filesep 'masks'];
refinedMasksOldDir = [analysisDir filesep 'refined_masks_OLD'];

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
refinedMasksDirFull = [refinedMasksOldDir filesep 'refined_masks_for_channel_1'];
refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tif']);
if isempty(refinedMaskFileListing)
    refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tiff']);
end

%get number of files
frames2rescue = frames2rescue(:)';
numFiles = length(frames2rescue);

%% Mask refinement

%refine masks using steerable line filter
particleInfo = [];
prctileUsed = NaN(numFiles,1);
for iFile = frames2rescue
    
    %read image and mask
    image = double(imread(fullfile(imageDir,imageFileListing(iFile).name)));
    mask0 = double(imread(fullfile(refinedMasksDirFull,refinedMaskFileListing(iFile).name)));
    
    %get relevant particle information
    if ~isempty(movieInfo)
        particleIndx = (iFile-1)*numSPTFrames;
        if iFile ~= numFiles
            particleInfo = movieInfo(particleIndx+1:particleIndx+numSPTFrames);
        end
    end
    
    %call refinement function
    [mask,prctileUsed(iFile)] = refineEdgeWithSteerableFilterSeed(image,mask0,...
        threshParam,gapCloseParam,doPlot,particleInfo,meanBkg);
    
    %store new mask
    imwrite(mask,fullfile(masksDirFull,maskFileListing(iFile).name),'tif');
    
    if prctileUsed(iFile)==-1
        disp(['bad mask Image ' num2str(iFile)]);
    end
    
end

%% ~~~ the end ~~~

