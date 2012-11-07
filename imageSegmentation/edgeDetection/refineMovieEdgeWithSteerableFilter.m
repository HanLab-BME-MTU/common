function prctileUsed = refineMovieEdgeWithSteerableFilter(MD,threshParam,...
    gapCloseParam,doPlot,meanBkg,channel2Dir)
%REFINEMOVIEEDGEWITHSTEERABLEFILTER replaces simple masks with masks refined using steerable line filtering
%
%SYNOPSIS prctileUsed = refineMovieEdgeWithSteerableFilter(MD,threshParam,...
%    gapCloseParam,doPlot,meanBkg,channel2Dir)
%
%INPUT  MD    : The movieData object as output by the cell masking and
%               windowing software. Before calling this code,
%               thresholdMovie and refineMovieMask must have been
%               processed. Otherwise the code will crash.
%       threshParam  : Structure with parameters for gradient thresholding:
%           .filterSigma    : Standard deviation for steerable filtering.
%                             Optional. Default: 1.5.
%           .prctile        : Percentile for thresholding.
%                             Optional. Default: [95 90 85 80].
%           .bandHalfWidth  : Half-width of band around original edge to
%                             look into for edge refinement.
%                             Input -1 to use whole image instead of a band.
%                             Optional. Default: 50 pixels.
%       gapCloseParam: Structure with parameters for edge gap closing:
%           .maxEdgePairDist: Maximum distance between edge segment pairs.
%                             Optional. Default: 5 pixels.
%           .factorContr    : Contribution of each factor to the edge gap
%                             closing cost. 4 entries for the factors:
%                             (1) distance,
%                             (2) angle between gradients,
%                             (3) angle between gradient and perpendicular
%                                 to centroid-centroid distance,
%                             (4) "edginess" score,
%                             Optional. Default: ones(1,4).
%       doPlot: 1 to plot masks in the end, 2 to also show edge progress,
%               0 to plot nothing. In final plot, refined masks shown in
%               green, original masks shown in blue. Note that this
%               refinement comes on top of the "refineMovieMask"
%               refinement.
%               Optional. Default: 0.
%       meanBkg  : Mean background intensity close to the cell edge.
%                  Optional. If not input, information not used.
%       channel2Dir: Directory where Channel 2 images, or some derivative
%                  thereof, are saved to use in edge refinement.
%                  Optional. If not input, information not used.
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
    threshParam.prctile = [95 90 85 80];
    threshParam.bandHalfWidth = 50;
else
    if ~isfield(threshParam,'fiterSigma')
        threshParam.filterSigma = 1.5;
    end
    if ~isfield(threshParam,'prctile')
        threshParam.prctile = [95 90 85 80];
    end
    if ~isfield(threshParam,'bandHalfWidth')
        threshParam.bandHalfWidth = 50;
    end
end

%get edge gap closing parameters
if nargin < 3 || isempty(gapCloseParam)
    gapCloseParam.maxEdgePairDist = 5;
    gapCloseParam.factorContr = ones(1,4);
else
    if ~isfield(gapCloseParam,'maxEdgePairDist')
        gapCloseParam.maxEdgePairDist = 5;
    end
    if ~isfield(gapCloseParam,'factorContr')
        gapCloseParam.factorContr = ones(1,4);
    end
end

%check whether/what to plot
if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end

%get background information
if nargin < 5 || isempty(meanBkg)
    meanBkg = [];
end

%check particle information
if nargin < 6 || isempty(channel2Dir)
    channel2Dir = [];
end

%get image and analysis directories
imageDir = MD.channels_.channelPath_;
analysisDir = [MD.movieDataPath_ filesep 'SegmentationPackage'];

%make new directories and copy old masks to them
masksDir = [analysisDir filesep 'masks'];
masksOldDir = [analysisDir filesep 'masks_OLD'];
mkdir(masksOldDir);
copyfile([masksDir filesep '*'],masksOldDir,'f');
refinedMasksDir = [analysisDir filesep 'refined_masks'];
refinedMasksOldDir = [analysisDir filesep 'refined_masks_OLD'];
mkdir(refinedMasksOldDir);
copyfile([refinedMasksDir filesep '*'],refinedMasksOldDir,'f');

%get image and mask file listings
imageFileListing = dir([imageDir filesep '*.tif']);
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.tiff']);
end
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.TIF']);
end
masksDirFull = [masksDir filesep 'masks_for_channel_1'];
maskFileListing = dir([masksDirFull filesep '*.tif']);
if isempty(maskFileListing)
    maskFileListing = dir([masksDirFull filesep '*.tiff']);
end
if isempty(maskFileListing)
    maskFileListing = dir([masksDirFull filesep '*.TIF']);
end
refinedMasksDirFull = [refinedMasksDir filesep 'refined_masks_for_channel_1'];
refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tif']);
if isempty(refinedMaskFileListing)
    refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tiff']);
end
if isempty(refinedMaskFileListing)
    refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.TIF']);
end

%get channel 2 file listings if relevant
if ~isempty(channel2Dir)
    channel2FileListing = dir([channel2Dir filesep '*.tif']);
    if isempty(channel2FileListing)
        channel2FileListing = dir([channel2Dir filesep '*.tiff']);
    end
    if isempty(channel2FileListing)
        channel2FileListing = dir([channel2Dir filesep '*.TIF']);
    end
end

%get number of files
numFiles = length(imageFileListing);

%% Mask refinement

wtBar = waitbar(0,'Please wait, refining edge with steerable filter ...'); 

%refine masks using steerable line filter
channel2Info = [];
prctileUsed = NaN(numFiles,1);
for iFile = 1 : numFiles
    
    waitbar(iFile/numFiles,wtBar,'Please wait, refining edge with steerable filter ...'); 
    
    %read image and mask
    image = double(imread(fullfile(imageDir,imageFileListing(iFile).name)));
    mask0 = double(imread(fullfile(refinedMasksDirFull,refinedMaskFileListing(iFile).name)));
    
    %get relevant channel2 information
    if ~isempty(channel2Dir)
        channel2Info = double(imread(fullfile(channel2Dir,channel2FileListing(iFile).name)));
    end
    
    %call refinement function
    [mask,prctileUsed(iFile)] = refineEdgeWithSteerableFilterSeed(image,mask0,...
        threshParam,gapCloseParam,doPlot,meanBkg,channel2Info);
    
    %store new mask
    imwrite(mask,fullfile(masksDirFull,maskFileListing(iFile).name),'tif');
    
    if prctileUsed(iFile)==-1
        disp(['bad mask Image ' num2str(iFile)]);
    end
    
end

close(wtBar)

%% ~~~ the end ~~~

