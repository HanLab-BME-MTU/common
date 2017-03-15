function detectMoviePointSources3D(movieDataOrProcess,varargin)
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

% Joy Xu / Sebastien Besson, July 2014

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('UseIntersection',false, @islogical);
% ip.addParamValue('EstimateSigma', false, @islogical); %%TODO
ip.parse(movieDataOrProcess,varargin{:});
ip.KeepUnmatched = true;
paramsIn = ip.Results.paramsIn;

% Get MovieData object and Process
[movieData, pointSourceDetProc3D] = getOwnerAndProcess(movieDataOrProcess,'PointSourceDetectionProcess3D',true);

assert(movieData.is3D, 'detectMoviePointSources3D is specifically designed for 3D images. Please use detectMoviePointSources for 2D images.');

%Parse input, store in parameter structure
p = parseProcessParams(pointSourceDetProc3D, paramsIn);


%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',pointSourceDetProc3D.getName());
else
    wtBar=-1;
end

% Reading various constants
nFrames = movieData.nFrames_;
nChanDet = numel(p.ChannelIndex);
imSize = movieData.imSize_;

%Find the  the segmentation process.
if isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex);
    p.MaskProcessIndex =movieData.getProcessIndex('MaskProcess',1,1);
end    

if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
    maskProc = movieData.processes_{p.MaskProcessIndex};
    if ~all(maskProc.checkChannelOutput(p.MaskChannelIndex))
        error('All channels must be segmented!')
    end
    
    %Create mask directory if several masks need to be merged
    if length(p.MaskChannelIndex) > 1 %&& p.UseIntersection
        close(wtBar); error('MaskIntersectionProcess not yet supported in 3D. Please sepecify a single MaskChannelIndex.'); 
%         %Get the indices of any previous mask intersection process
%         iMaskIntProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
%         
%         %If the process doesn't exist, create it
%         if isempty(iMaskIntProc)
%             iMaskIntProc = numel(movieData.processes_)+1;
%             movieData.addProcess(MaskIntersectionProcess(movieData,p.OutputDirectory));
%         end
%         maskIntProc = movieData.processes_{iMaskIntProc};
%         
%         %Set up the parameters for mask intersection
%         maskIntParams.ChannelIndex = p.MaskChannelIndex;
%         maskIntParams.SegProcessIndex = p.MaskProcessIndex;
%         
%         parseProcessParams(maskIntProc,maskIntParams);
%         maskIntProc.run;
%         
%         %Then use this mask process
%         maskProc = maskIntProc;
            
    end 
    
    % Get mask directory and names
    maskDir = maskProc.outFilePaths_(p.MaskChannelIndex);
    
end

%Check the input processes if any and get loader handles for each channel
imDirs = cell(1,nChanDet);
imLoader = cell(1,nChanDet);
for j = 1:nChanDet
    
    if p.InputImageProcessIndex(j) > 0
        %Check the specified input process
        assert(movieData.processes_{p.InputImageProcessIndex(j)}.checkChannelOutput(p.ChannelIndex(j)),['No valid output for input process specified for channel ' num2str(p.ChannelIndex(j))]);
        imDirs{p.ChannelIndex(j)} = movieData.processes_{p.InputImageProcessIndex(j)}.outFilePaths_{p.ChannelIndex(j)};
        imLoader{p.ChannelIndex(j)} = @(f)(movieData.processes_{p.InputImageProcessIndex(j)}.loadChannelOutput(p.ChannelIndex(j),f));                    
        
    else
        %If raw data specified
        imDirs{p.ChannelIndex(j)} = movieData.channels_(p.ChannelIndex(j)).channelPath_;
        imLoader{p.ChannelIndex(j)} = @(f)(movieData.channels_(p.ChannelIndex(j)).loadStack(f));
        
    end
    
end



% Set up the input directories
inFilePaths = cell(1,numel(movieData.channels_));
for j = 1:numel(p.ChannelIndex)
    inFilePaths{1,p.ChannelIndex(j)} = imDirs{p.ChannelIndex(j)};
    if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
        inFilePaths{2,p.ChannelIndex(j)} = maskDir;
    end
end
pointSourceDetProc3D.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory)
pointSourceDetProc3D.setOutFilePaths(outFilePaths);

%Get ROI mask if any.
%roiMask = movieData.getROIMask;

%% --------------- Add optional auio-estimation of PSF sigma ---------------%%% 



%% --------------- Point source detection ---------------%%% 

disp('Starting detecting diffraction-limited objects');

logMsg = @(chan) ['Please wait, detecting diffraction-limited objects for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;
nj = 0;
for i = 1:numel(p.ChannelIndex)
    
    iChan = p.ChannelIndex(i);
    % Log display
    if ishandle(wtBar), waitbar(nj/nTot,wtBar,sprintf(logMsg(iChan)));  end
    disp(logMsg(iChan))
    disp(imDirs{1,iChan});
    if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
        disp(sprintf('Using mask from: %s', maskDir{1}));
    end
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    %Set up parameter structure for detection on this channel
    detP = splitPerChannelParams(p, iChan);
    
    for j = 1:nFrames
                
        % loading the entire stack
        vol = imLoader{p.ChannelIndex(i)}(j);
        
        %!!!!! fix the mask for 3d data!!!!!%
        if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
            currMask = maskProc.loadChannelOutput(p.MaskChannelIndex,j); %& roiMask(:,:,j); %istack is the stack index!!!! test for validity!!!!
            detP.Mask =  currMask;
        else
            detP.Mask = true([movieData.imSize_(1), movieData.imSize_(2), movieData.zSize_]);%roiMask(:,:,j);
        end    
        
        % Call main detection function       
        pstruct = pointSourceDetection3D(vol,p.filterSigma(:,iChan),detP);


        % add xCoord, yCoord, amp fields for compatibilty  with tracker
        if ~isempty(pstruct)
            pstruct.xCoord = [pstruct.x' pstruct.x_pstd'];
            pstruct.yCoord = [pstruct.y' pstruct.y_pstd'];
            pstruct.zCoord = [pstruct.z' pstruct.z_pstd'];
            pstruct.amp = [pstruct.A' pstruct.A_pstd'];        
            pstruct.sigmaX = [pstruct.s' pstruct.s_pstd'];
            pstruct.sigmaY = [pstruct.s' pstruct.s_pstd'];
            pstruct.sigmaZ = [pstruct.s' pstruct.s_pstd'];
            pstruct.bkg = [pstruct.c' pstruct.c_pstd'];
        end
        
        if ~isempty(pstruct) && ~exist('allFields','var')
            %Initialize the structure
            allFields = fieldnames(pstruct);
            allFields = vertcat(allFields', arrayfun(@(x)([]),1:numel(allFields),'Unif',false));            
        end
        if j == 1 && exist('allFields','var')
            movieInfo(1:nFrames)= struct(allFields{:});                        
        end
        
        if ~isempty(pstruct)
            movieInfo(j) = pstruct; %#ok<AGROW>
        end
        
        if ishandle(wtBar) && (nFrames <  20 || mod(j,5)==1 )
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    if ~exist('movieInfo','var')
        %in the case that no channels/frames had detected points
        movieInfo = [];
    end
    save(outFilePaths{1,iChan}, 'movieInfo'); 
    clear movieInfo;
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished detecting diffraction-limited objects!')

movieData.save;


