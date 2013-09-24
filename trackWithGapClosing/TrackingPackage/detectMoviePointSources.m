function detectMoviePointSources(movieData,varargin)
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
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('PointSourceDetectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(PointSourceDetectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
pointSourceDetProc = movieData.processes_{iProc};
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
bitDepth = movieData.camBitdepth_;
nFrames = movieData.nFrames_;
maxIntensity =(2^bitDepth-1);

%Find the  the segmentation process
if isempty(p.MaskProcessIndex) 
    p.MaskProcessIndex =movieData.getProcessIndex('MaskProcess',1,1);
end    

if isempty(p.MaskChannelIndex)
    p.MaskChannelIndex = p.ChannelIndex;
end

if ~isempty(p.MaskProcessIndex)
    maskProc = movieData.processes_{p.MaskProcessIndex};
    if ~all(maskProc.checkChannelOutput(p.MaskChannelIndex))
        error('All channels must be segmented!')
    end
    
    %Create mask directory if several masks need to be merged
    if length(p.MaskChannelIndex) >1
        %Get the indices of any previous mask intersection process
        iMaskIntProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
        
        %If the process doesn't exist, create it
        if isempty(iMaskIntProc)
            iMaskIntProc = numel(movieData.processes_)+1;
            movieData.addProcess(MaskIntersectionProcess(movieData,p.OutputDirectory));
        end
        maskIntProc = movieData.processes_{iMaskIntProc};
        
        %Set up the parameters for mask transformation
        maskIntParams.ChannelIndex = p.MaskChannelIndex;
        maskIntParams.SegProcessIndex = p.MaskProcessIndex;
        
        parseProcessParams(maskIntProc,maskIntParams);
        maskIntProc.run;
        
        % Get mask directory and names
        maskDir = maskIntProc.outFilePaths_{1};
    else
        maskDir = maskProc.outFilePaths_{p.MaskChannelIndex};
    end    
end

% Set up the input directories
inFilePaths = cell(1,numel(movieData.channels_));
for j = p.ChannelIndex
    inFilePaths{1,j} = imDirs{j};
    if ~isempty(p.MaskProcessIndex)
        inFilePaths{2,j} = maskDir;
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
    if ~isempty(p.MaskProcessIndex)
        disp(sprintf('Using mask from: %s', maskDir));
    end
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    for j= 1:nFrames

        currImage = double(movieData.channels_(iChan).loadImage(j))/maxIntensity; 
        if ~isempty(p.MaskProcessIndex)
            currMask = maskProc.loadChannelOutput(p.MaskChannelIndex(1),j);
            p.Mask =  currMask;            
        else
            p.Mask = [];
        end    
        
        % Call main detection function
        pstruct = pointSourceDetection(currImage,p.filterSigma,p);


        % add xCoord, yCoord, amp fields for compatibilty  with tracker
        pstruct.xCoord = [pstruct.x' pstruct.x_pstd'];
        pstruct.yCoord = [pstruct.y' pstruct.y_pstd'];
        pstruct.amp = [pstruct.A' pstruct.A_pstd'];
        
        if j == 1
            movieInfo(1:nFrames) = pstruct;
        else
            movieInfo(j) = pstruct;
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

