function computeMovieMIP(MD, varargin)
%computeMovieMIP MovieData-Process framework wrapped printMIP function
% Andrew R. Jamieson, Aug. 2017
% Also See printMIP                                    
    
ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(MD, varargin{:});
paramsIn = ip.Results.paramsIn;

if(MD.zSize_==1)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

% Gather process and MD info/parameters
[movieData1, thisProcess, iProc] = getOwnerAndProcess(MD, 'ComputeMIPProcess', true);
assert(MD == movieData1);

%Parse input, store in parameter structure
p = parseProcessParams(thisProcess, paramsIn);

% ============= Configure InputPaths. ================
inFilePaths = cell(1, numel(MD.channels_));
for j = p.ChannelIndex
    inFilePaths{1,j} = MD.getChannelPaths{j};
end
thisProcess.setInFilePaths(inFilePaths);

% ================[OUTPUT]===========================
mkClrDir(p.OutputDirectory);
outFilePaths = cell(3, numel(MD.channels_));
for i = p.ChannelIndex;    
    outFilePaths{1,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'XY'];
    outFilePaths{2,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'ZY'];
    outFilePaths{3,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'ZX'];
    outFilePaths{4,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'three'];
    outFilePaths{5,i} = [p.OutputDirectory filesep 'ch' num2str(i)];
    mkClrDir(outFilePaths{1,i});
    mkClrDir(outFilePaths{2,i});
    mkClrDir(outFilePaths{3,i});
    mkClrDir(outFilePaths{4,i});
end
thisProcess.setOutFilePaths(outFilePaths);

ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;
minIntensityNorm = [];
maxIntensityNorm = [];

% prep for parfor linearization 
pRunSet = {}; cell(MD.nFrames_ * length(p.ChannelIndex));
i = 1;
for chIdx = p.ChannelIndex
    for frameIdx = 1:MD.nFrames_
        pRunSet{i}.chIdx = chIdx; %struct('chIdx', chIdx)
        pRunSet{i}.frameIdx = frameIdx; %struct('chIdx', chIdx)
        i = i + 1;
    end
end

parfor pfi = 1:length(pRunSet)  
% for chIdx = p.ChannelIndex
    chIdx = pRunSet{pfi}.chIdx;
    frameIdx = pRunSet{pfi}.frameIdx;
    
    savePathXY = outFilePaths{1, chIdx};
    savePathZY = outFilePaths{2, chIdx};
    savePathZX = outFilePaths{3, chIdx};
    savePathThree = outFilePaths{4, chIdx};
    savePath = outFilePaths{5, chIdx};
    if ~isdir(savePath) || ~isdir(savePathXY) || ~isdir(savePathZY) || ~isdir(savePathZX) || ~isdir([savePath filesep 'Three'])
        mkdirRobust([savePath]);
        mkdirRobust([savePathXY]);
        mkdirRobust([savePathZY]);
        mkdirRobust([savePathZX]);
        mkdirRobust([savePathThree]);
    end

    XYFilesPattern = [savePathXY filesep 'XY_frame_nb%04d.png'];
    YZFilesPattern = [savePathZY filesep 'ZY_frame_nb%04d.png'];
    XZFilesPattern = [savePathZX filesep 'ZX_frame_nb%04d.png'];
    ThreeFilesPattern = [savePathThree filesep 'Three_frame_nb%04d.png'];

    vol = MD.getChannel(chIdx).loadStack(1);
    minIntensityNorm = min(vol(:));
    maxIntensityNorm = max(vol(:));

%     parfor frameIdx = 1:MD.nFrames_
        
        vol = MD.getChannel(chIdx).loadStack(frameIdx);  
        [maxXY, maxZY, maxZX, three] = computeMIPs(vol, ZXRatio, minIntensityNorm, maxIntensityNorm);
        
        % save the maximum intensity projections
        imwrite(maxXY, sprintfPath(XYFilesPattern, frameIdx), 'Compression', 'none');
        imwrite(maxZY, sprintfPath(YZFilesPattern, frameIdx), 'Compression', 'none');
        imwrite(maxZX, sprintfPath(XZFilesPattern, frameIdx), 'Compression', 'none');
        imwrite(three, sprintfPath(ThreeFilesPattern, frameIdx), 'Compression', 'none');
%     end
end

parfor chIdx = p.ChannelIndex
    % savePath = outFilePaths{1,chIdx};
    ThreeFilesPattern = [outFilePaths{4, chIdx} filesep 'Three_frame_nb%04d.png'];
    threeVideo = VideoWriter([outFilePaths{5, chIdx} filesep 'threeMontage.avi']);
    % myVideo.FrameRate = 4;  % Default 30
    % myVideo.Quality = 90;    % Default 75

    open(threeVideo);
    for frameIdx = 1:MD.nFrames_
        three=imread(sprintfPath(ThreeFilesPattern, frameIdx));
        writeVideo(threeVideo,three)
    end
    close(threeVideo);
end