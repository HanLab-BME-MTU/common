function detectMoviePointSources(movieDataOrProcess,varargin)
% detectMoviePointSource detect diffraction-limited objects in a movie
%
% detectMoviePointSources
%
% SYNOPSIS detectMoviePointSources(movieData,paramsIn)
%
% INPUT
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT

% Sebastien Besson, Sep 2011 (last modified Mar 2013)

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('UseIntersection',true,@islogical);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get MovieData object and Process
[movieData, pointSourceDetProc] = getOwnerAndProcess(movieDataOrProcess,'PointSourceDetectionProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(pointSourceDetProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',pointSourceDetProc.getName());
else
    wtBar=-1;
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
nFrames = movieData.nFrames_;

%Find the  the segmentation process.
if isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex);
    p.MaskProcessIndex =movieData.getProcessIndex('MaskProcess',1,1);
end

if ~isempty(p.MaskChannelIndex) && numel(p.MaskChannelIndex) ~= numel(p.ChannelIndex)
    error('If masks are applied you must specify one mask channel per detection channel!')
end

if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
    maskProc = movieData.processes_{p.MaskProcessIndex};
    if ~all(maskProc.checkChannelOutput(p.MaskChannelIndex))
        error('All channels must be segmented!')
    end
    
    %Create mask directory if several masks need to be merged
    if length(p.MaskChannelIndex) > 1 && p.UseIntersection
        %Get the indices of any previous mask intersection process
        iMaskIntProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
        
        %If the process doesn't exist, create it
        if isempty(iMaskIntProc)
            iMaskIntProc = numel(movieData.processes_)+1;
            movieData.addProcess(MaskIntersectionProcess(movieData,p.OutputDirectory));
        end
        maskIntProc = movieData.processes_{iMaskIntProc};
        
        %Set up the parameters for mask intersection
        maskIntParams.ChannelIndex = p.MaskChannelIndex;
        maskIntParams.SegProcessIndex = p.MaskProcessIndex;
        
        parseProcessParams(maskIntProc,maskIntParams);
        maskIntProc.run;
        
        %Then use this mask process
        maskProc = maskIntProc;
        
    end
    
    % Get mask directory and names
    maskDir = maskProc.outFilePaths_(p.MaskChannelIndex);
    
end

% Set up the input directories
inFilePaths = cell(1,numel(movieData.channels_));
for j = 1:numel(p.ChannelIndex)
    inFilePaths{1,p.ChannelIndex(j)} = imDirs{p.ChannelIndex(j)};
    if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
        inFilePaths{2,p.ChannelIndex(j)} = maskDir{j};
    end
end
pointSourceDetProc.setInFilePaths(inFilePaths);

% Set up the output directories
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory)
pointSourceDetProc.setOutFilePaths(outFilePaths);

% Get ROI mask if any.
roiMask = movieData.getROIMask;

%% --------------- Point source detection ---------------%%%

disp('Starting detecting diffraction-limited objects');

logMsg = @(chan) ['Please wait, detecting diffraction-limited objects for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;

for i = 1:numel(p.ChannelIndex)
    
    iChan = p.ChannelIndex(i);
    % Log display
    if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg(iChan)));  end
    disp(logMsg(iChan))
    disp(imDirs{1,iChan});
    if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
        fprintf(1, 'Using mask from: %s', maskDir{i});
    end
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    %Set up parameter structure for detection on this channel
    detP = splitPerChannelParams(p, iChan);
    detP.Mask = roiMask(:,:,1);
    
    % Initialize a movieInfo structure with the minimum amount of fields
    clear movieInfo
    movieInfo(1 : nFrames)= struct('xCoord', [], 'yCoord', [],...
        'amp', [], 'sigmaX', [], 'sigmaY', [], 'bkg', []);
    
    parfor j= 1:nFrames
%     for j= 1:nFrames
        
        currImage = double(movieData.channels_(iChan).loadImage(j));
        if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
            currMask = maskProc.loadChannelOutput(p.MaskChannelIndex(i),j) & roiMask(:,:,j);
%             detP.Mask =  currMask;
        else
            currMask = roiMask(:,:,j);
        end
%         detP.Mask = currMask;
        
        % Call main detection function
        pstruct{j} = pointSourceDetection(currImage,p.filterSigma(iChan),detP,'Mask',currMask);
    end
    for j= 1:nFrames
        if ~isempty(pstruct{j}) &&...
                ~isequal(fieldnames(pstruct{j}), fieldnames(movieInfo(j)))
            allFields = fieldnames(pstruct{j});
            for iField = 1:numel(allFields)
                movieInfo(j).(allFields{iField}) = pstruct{j}.(allFields{iField}); %#ok<AGROW>
            end
            movieInfo = orderfields(movieInfo);
        end
        % add xCoord, yCoord, amp fields for compatibilty  with tracker
        if ~isempty(pstruct{j})
            
            pstructOrdered = pstruct{j};
            pstructOrdered.xCoord = [pstruct{j}.x' pstruct{j}.x_pstd'];
            pstructOrdered.yCoord = [pstruct{j}.y' pstruct{j}.y_pstd'];
            pstructOrdered.amp = [pstruct{j}.A' pstruct{j}.A_pstd'];
            pstructOrdered.sigmaX = [pstruct{j}.s' pstruct{j}.s_pstd'];
            pstructOrdered.sigmaY = [pstruct{j}.s' pstruct{j}.s_pstd'];
            pstructOrdered.bkg = [pstruct{j}.c' pstruct{j}.c_pstd'];
            
            movieInfo(j) = orderfields(pstructOrdered); %#ok<AGROW>
        end
        
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    save(outFilePaths{1,iChan}, 'movieInfo');
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished detecting diffraction-limited objects!')

