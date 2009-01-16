function [image,reader,metadata] = imreadLoci(fileName,tRange,cRange,zRange,verbose)
%IMREADLOCI reads images using the bio-formats java library
%
%SYNOPSIS: [image,reader,metadata] = imreadLoci(fileName,tRange,cRange,zRange,verbose)
%
% INPUT fileName: Full filename including path of the image file you try to
%           open. If empty, imreadLoci will open a dialog.
%           Alternatively, you can pass a reader object
%       tRange/cRange/zRange: vectors with timepoint(s), color(s) and
%           z-slice(s) to read. If tRange is -1, only metadata
%           is read. Default: read all.
%       verbose : if 1, progressText is shown. Default: 1.
%
% OUTPUT image : image array
%        reader : reader object. Pass this back to the function to speed up
%                 loading.
%        metadata : struct with the following currently supported fields
%           imageSize : imageSize in voxels. Ordering XYZTC
%           microscopySettings : base microscopy settings that are required
%               for quantitative image analysis. Struct with fields
%               .voxelSize : voxelSize XYZ in microns
%               .timeLapse : average (robust) time between time points
%               .NA        : numerical aperture of the objective
%               .wavelength: emission wavelength in microns
%           microscopySettingsExt : extended microscopy settings. The
%               following fields are currently filled in automatically
%               .timestamps     : individual timestamps of z-slices
%               .timeLapseStats : mean,sd,%outliers of time between frames
%           experimentInfo : information about the experiment. The
%               following fields are currently being filled in
%               automatically if possible
%               .date       : creation date
%               .creator    : creator of movie
%               .waveLabels : channel labels
%
% REMARKS:  This function is not quite finished yet.
%           You need loci_tools.jar on your java path.
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-Sep-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign empty to allow return below
image = [];
reader = [];
metadata = [];

%=========================
% TEST INPUT
%=========================

% check for fileName and load if necessary
if nargin == 0 || isempty(fileName)
    [fileName,pathName] = uigetfile({'*.*','All Files'},'Please select an image file');
    if fileName == 0
        % user aborted. Exit gracefully
        
        disp('No file selected')
        return
    else
        fileName = fullfile(pathName,fileName);
    end
else
    % check if file exists
    if isa(fileName,'loci.formats.ChannelSeparator')
        reader = fileName;
        fileName = reader.getCurrentFile.toCharArray';
    elseif ~exist(fileName,'file')
        error('file %s not found',fileName);
    end
end

% check for ranges
if nargin < 2 || isempty(tRange)
    tRange = 0;
end
if nargin < 3 || isempty(cRange)
    cRange = 0;
end
if nargin < 4 || isempty(zRange)
    zRange = 0;
end
if nargin < 5 || isempty(verbose)
    verbose = true;
end

% check whether to read movie
readMovie = ~(tRange == -1);

% check whether to read metadata
readMetadata = nargout > 2;

% check whether reader is already defined
defineReader = isempty(reader);


%============================



%============================
% CREATE OBJECTS
%============================

% find image name
[imagePath, imageName] = fileparts(fileName);

if defineReader
    
    % create reader object
    reader = loci.formats.ChannelSeparator;
    
    % allow reading metadata
    if readMetadata
        store = loci.formats.MetadataTools.createOMEXMLMetadata;
        
        % add metaData-container to reader
        reader.setMetadataStore(store);
    end
    
    % supply the image
    %progressText(0,sprintf('reading properties for %s',imageName))
    reader.setId(fileName);
    
end

%======================
%% READ METADATA
%======================

% read image size - this we can get from reader directly
imSize = zeros(1,5);
imSize(1) = reader.getSizeX;
imSize(2) = reader.getSizeY;
imSize(3) = reader.getSizeZ;
imSize(4) = reader.getSizeT;
imSize(5) = reader.getSizeC;

% read metadata
if readMetadata
    
    if verbose
        progressText(0,sprintf('reading metadata for %s',imageName))
    end
    
    % get the metadata object
    meta = reader.getMetadataStore;
    
    metadata.imageSize = imSize;
    
    % read base microscopy data
    metadata.microscopySettings = struct('voxelSize',NaN(1,3),'timeLapse',NaN,...
        'NA',NaN,'wavelength',NaN(1,5));
    % also set up extended microscopy settings, as well as experimenter
    % info. Not all of the fields are currently being read
    metadata.microscopySettingsExt = struct('exposureTime',NaN(1,5),...
        'microscopeType','','magnification',NaN,'auxMag',NaN,...
        'binning',NaN,'gain',NaN,'camReadSpeed',NaN,...
        'aperture',{{NaN,NaN}},'ndFilter',NaN(1,5),...
        'laserPower',NaN(1,5),'timestamps',[],'timeLapseStats',NaN(1,3));
    metadata.experimentInfo = struct('date',[],'creator',[],'waveLabels',{cell(1,5)});
    
    % read voxelSize.
    dim = 'XYZ';
    for d = 1:3
        tmp = eval(sprintf('meta.getDimensionsPhysicalSize%s(0,0);',dim(d)));
        if ~isempty(tmp)
            metadata.microscopySettings.voxelSize(d) = tmp;
        end
    end
    % read timing
    % Unfortunately, there is no field like timeLapse in the metaObj.
    % Read individual timestamps and estimate tL from there. However, this
    % is somewhat tricky because of dimension order in which the data is
    % stored
    timeStamp = NaN(prod(imSize(3:5)),1);
    for t=1:length(timeStamp)
        tmp = meta.getPlaneTimingDeltaT(0,0,t-1);
        if ~isempty(tmp)
            timeStamp(t) = tmp.doubleValue;
        end
    end
    
    if any(isnan(timeStamp)) || all(timeStamp == 0)
        warning('IMREADLOCI:CANNOTREADMETADATA',...
            'timestamps could not be read')
    else
        
        dimOrder = reader.getDimensionOrder.toCharArray';
        if strmatch(dimOrder(end),'T')
            % time is last in dimension order.
            tmp = reshape(timeStamp,imSize(3)*imSize(5),[]);
            tmp = diff(tmp,1,2);
            % get dt using robustMean
            [dt,rs,ii] = robustMean(tmp(:));
            op = 100*(1-length(ii)/numel(tmp));
            metadata.microscopySettings.timeLapse = dt;
            metadata.microscopySettingsExt.timeLapseStats = [dt,rs,op];
        else
            warning('IMREADLOCI:CANNOTREADMETADATA',...
                'dimension order %s is not supported yet for reading timestamp',dimOrder);
        end
    end
    % read NA
    tmp = meta.getObjectiveLensNA(0,0);
    if ~isempty(tmp)
        metadata.microscopySettings.NA = tmp.doubleValue;
    end
    
    % read wavelength % store labels
    for w = 1:imSize(5)
        tmp = meta.getLogicalChannelEmWave(0,w-1);
        if ~isempty(tmp)
            % convert to microns!
            metadata.microscopySettings.wavelength(w) = tmp.doubleValue/1000;
        end
        tmp = meta.getLogicalChannelName(0,w-1);
        if ~isempty(tmp)
            metadata.experimentInfo.waveLabels{w} = tmp.toCharArray';
        end
    end
    
    % -------- read extended info here -------
    
    
    % read experiment info
    d = dir(fileName);
    metadata.experimentInfo.date = d.date;
    % for now, only store one experimentator
    fn = meta.getExperimenterFirstName(0);
    ln = meta.getExperimenterLastName(0);
    em = meta.getExperimenterEmail(0);
    if ~isempty(fn)
        fn = fn.toCharArray';
    end
    if ~isempty(ln)
        ln = ln.toCharArray';
    end
    if ~isempty(em)
        em = em.toCharArray';
    end
    if ~(isempty(fn) && isempty(ln) && isempty(em))
        metadata.experimentInfo.creator = sprintf('%s %s %s',fn,ln,em);
    end
    
    if verbose
        progressText(1)
    end
    
end

%========================
%% READ MOVIE
%========================

if readMovie
    
    % check for size and range
    if isscalar(zRange) && zRange == 0
        zRange = 1:imSize(3);
    end
    if isscalar(tRange) && tRange == 0
        tRange = 1:imSize(4);
    end
    if isscalar(cRange) && cRange == 0
        cRange = 1:imSize(5);
    end
    
    % check for out of range
    if any(zRange>imSize(3))
        error('there are only %i z-slices',imSize(3))
    end
    if any(tRange>imSize(4))
        error('there are only %i timepoints',imSize(4))
    end
    if any(cRange>imSize(5))
        error('there are only %i channels',imSize(5))
    end
    
    nZ = length(zRange);
    nT = length(tRange);
    nC = length(cRange);
    
    % preassign arrays
    image = zeros(imSize(1),imSize(2),nZ,nT,nC);
    imagePlane = zeros(imSize(1:2));
    
    
    % read images and store
    if verbose
        progressText(0,sprintf('loading %s',imageName));
        nImages = nZ*nT*nC;
        imageCt = 0;
    end
    
    
    for c = cRange(:)'
        for t = tRange(:)'
            for z = zRange(:)'
                % get current index
                idx = reader.getIndex(z-1,c-1,t-1);
                % read image
                img = reader.openImage(idx);
                % convert Java BufferedImage to MATLAB image
                imagePlane(:) = img.getData.getPixels(0, 0, imSize(1), imSize(2), []);
                             
                % store
                image(:,:,z,t,c) = imagePlane;
                
                % report progress
                if verbose
                    imageCt = imageCt + 1;
                    progressText(imageCt/nImages);
                end
            end
        end
    end
    
end
