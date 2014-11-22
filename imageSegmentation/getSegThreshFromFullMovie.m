function thresholdValue = getSegThreshFromFullMovie(MD,filterSigma,incAftMin,showPlot)

%% Input

%ask user for images
% if nargin < 1 || isempty(firstImageFile)
%     [fName,dirName] = uigetfile('*.tif','Specify first image in stack');
% else
%     if iscell(firstImageFile)
%         [fpath,fname,fno,fext]=getFilenameBody(firstImageFile{1});
%         dirName=[fpath,filesep];
%         fName=[fname,fno,fext];
%     elseif ischar(firstImageFile)
%         [fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
%         dirName=[fpath,filesep];
%         fName=[fname,fno,fext];
%     end
% end
% 
% %if input is valid ...
% if(isa(fName,'char') && isa(dirName,'char'))
%     
%     %get all file names in stack
%     outFileList = getFileStackNames([dirName,fName]);
%     numFrames = length(outFileList);
%     
% else %else, exit
%     
%     disp('--getSegThreshFromFullMovie: Bad file selection');
%     return
%     
% end

imageDir = MD.channels_.channelPath_;
analysisDir = MD.movieDataPath_;

imageFileListing = dir([imageDir filesep '*.tif']);
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.tiff']);
end
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.TIF']);
end
numFrames = length(imageFileListing);

if nargin < 2 || isempty(filterSigma)
    filterSigma = 0;
end

if nargin < 3 || isempty(incAftMin)
    incAftMin = 0;
end

if nargin < 4 || isempty(showPlot)
    showPlot = 0;
end

%% Threshold

%read all images
% imageStack = imread(outFileList{1});
imageStack = imread(fullfile(imageDir,imageFileListing(1).name));
imageStack = repmat(imageStack,[1 1 numFrames]);
for iFrame = 2 : numFrames
    %     imageStack(:,:,iFrame) = imread(outFileList{iFrame});
    imageStack(:,:,iFrame) = imread(fullfile(imageDir,imageFileListing(iFrame).name));
end
imageStack = double(imageStack);

%filter images
if filterSigma > 0
    for iFrame = 1 : numFrames
        imageStack(:,:,iFrame) = filterGauss2D(imageStack(:,:,iFrame),filterSigma);
    end
end

%START COPY-PASTE FROM THRESHOLDFLUORESCENCEIMAGE

%Get histogram, using Jonas' automatic bin number selection & smoothing
[~,~,histSpline] = optimalHistogram(imageStack(:),'smooth');

%Find the location of extrema in the histogram
histExtrema = fnzeros(fnder(histSpline));

%Get rid of the 'fake' extrema sometimes produced at beginning and end of
%distrubution by spline form.
histExtrema = histExtrema(1,:); %remove the intervals
histExtrema = histExtrema((histExtrema ~= ... %These will always be at first or last breaks in spline
            histSpline.breaks(1)) & (histExtrema ~= histSpline.breaks(end)));
histExtVals = fnval(histSpline,histExtrema);

%Determine whether each extrema is maximum or minimum
isMax = fnval(fnder(histSpline,2),histExtrema) < 0;

%Find the lowest-intensity maximum, assume this is the background peak.
iBackMax = find(isMax,1,'first');

%Find the first minimum after this maximum. This is used as the threshold.
iSep = iBackMax + 1;

thresholdValue = histExtrema(iSep);
minVal = histExtVals(iSep);

%END COPY-PASTE FROM THRESHOLDFLUORESCENCEIMAGE

%modify threshold by going away from the minimum
%this tries to handle flat minima
if incAftMin ~= 0
    minVal = minVal * (1 + incAftMin);
    threshTmp = thresholdValue:0.5:histExtrema(iSep+1);
    histVal = ppval(histSpline,threshTmp);
    histVal = abs(histVal - minVal);
    thresholdValue = threshTmp(histVal==min(histVal));
end

if showPlot
    figure, hold on
    optimalHistogram(imageStack(:))
    fnplt(histSpline,'r')        
    plot(histExtrema,histExtVals,'oc')
    plot(thresholdValue,minVal,'xg')
end

save(fullfile(analysisDir,'thresholdParamValue'),'filterSigma','incAftMin','thresholdValue');

%% ~~~ the end ~~~

