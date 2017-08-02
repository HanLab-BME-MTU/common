function printMIP(MD, varargin)
% Title says it all -- PR, augmented from Meghan D. 2015
% Input Parameter arguments:
%   Channel<numeric>: specifies which channel to MIP process
%                     default (non provided) will combine all into montage image.
%                                     
    
ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addOptional('Channel', [], @(x) isnumeric(x) && ismember(x, 1:length(MD.channels_)))
ip.parse(MD, varargin{:});

selectedChannel = ip.Results.Channel;
%progressText(0,'Print MIP');

% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

if(MD.zSize_==1)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

if isempty(selectedChannel)
    savePath=[MD.outputDirectory_ filesep 'MIP'];
else
    savePath=[MD.outputDirectory_ filesep 'MIP' filesep 'ch' num2str(selectedChannel)];
end

if ~isdir(savePath) || ~isdir([savePath filesep 'XY']) || ~isdir([savePath filesep 'ZY']) || ~isdir([savePath filesep 'ZX']) || ~isdir([savePath filesep 'Three'])
    system(['mkdir -p "' savePath '"']);
    system(['mkdir -p "' [savePath filesep 'XY'] '"' ]);
    system(['mkdir -p "' [savePath filesep 'ZY'] '"']);
    system(['mkdir -p "' [savePath filesep 'ZX'] '"']);
    system(['mkdir -p "' [savePath filesep 'Three'] '"']);
end

ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;
minIntensityNorm = [];
maxIntensityNorm = [];

if isempty(selectedChannel)
    for chIdx = 1:length(MD.channels_)
        vol = MD.getChannel(chIdx).loadStack(1);
        minIntensityNorm = [minIntensityNorm min(vol(:))];
        maxIntensityNorm = [maxIntensityNorm max(vol(:))];
    end
else
    vol = MD.getChannel(selectedChannel).loadStack(1);
    minIntensityNorm = min(vol(:));
    maxIntensityNorm = max(vol(:));
end


if isempty(selectedChannel)
    nameCells = MD.getChannel(1).getImageFileNames;
else
    nameCells = MD.getChannel(selectedChannel).getImageFileNames;
end

parfor frameIdx = 1:MD.nFrames_

    if isempty(selectedChannel)
        maxXY=[]; maxZY=[]; maxZX=[]; three=[];
        for chIdx = 1:length(MD.channels_)
            vol = MD.getChannel(chIdx).loadStack(frameIdx);  
            [cmaxXY, cmaxZY, cmaxZX, cthree] = computeMIPs(vol, ZXRatio, minIntensityNorm(chIdx), maxIntensityNorm(chIdx));
            maxXY = [maxXY cmaxXY];
            maxZY = [maxZY cmaxZY];
            maxZX = [maxZX cmaxZX];
            three = [three cthree];        
        end
    else
        chIdx = selectedChannel;
        vol = MD.getChannel(chIdx).loadStack(frameIdx);  
        [maxXY, maxZY, maxZX, three] = computeMIPs(vol, ZXRatio, minIntensityNorm, maxIntensityNorm);
    end
    
    % save the maximum intensity projections
    imwrite(maxXY, [savePath filesep 'XY' filesep 'XY_' nameCells{frameIdx} ], 'Compression', 'none');
    imwrite(maxZY, [savePath filesep 'ZY' filesep 'ZY_'  nameCells{frameIdx}], 'Compression', 'none');
    imwrite(maxZX, [savePath filesep 'ZX' filesep 'ZX_'  nameCells{frameIdx}], 'Compression', 'none');
    imwrite(three, [savePath filesep 'Three' filesep nameCells{frameIdx}], 'Compression', 'none');

end

threeVideo = VideoWriter([savePath filesep 'three.avi']);
myVideo.FrameRate = 4;  % Default 30
myVideo.Quality = 90;    % Default 75

open(threeVideo);
for frameIdx = 1:MD.nFrames_
    % save the maximum intensity projections
    three = imread([savePath filesep 'Three' filesep nameCells{frameIdx}]);
    writeVideo(threeVideo, three)
end
close(threeVideo);