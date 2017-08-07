function process=printMIP(MD)
% Title says it all -- PR, augmented from Meghan D. 2015

ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.parse(MD);

%progressText(0,'Print MIP');

% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

if(MD.zSize_==1)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

% processMIP=ExternalProcess(MD,'printMIP');
%     p.processSingleProj.setOutFilePaths({[outputDirSingleProj filesep 'XY' filesep 'frame_nb%04d.png'], ...
%         [outputDirSingleProj filesep 'YZ' filesep 'frame_nb%04d.png'], ...
%         [outputDirSingleProj filesep 'XZ' filesep 'frame_nb%04d.png'], ...
%         [outputDirSingleProj filesep 'limits.mat']});
%     frameNb=MD.nFrames_;
%     save([outputDirSingleProj filesep 'limits.mat'],'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');

savePath=[MD.outputDirectory_ filesep 'MIP'];
if ~isdir(savePath) || ~isdir([savePath filesep 'XY']) || ~isdir([savePath filesep 'ZY']) || ~isdir([savePath filesep 'ZX']) || ~isdir([savePath filesep 'three'])
    mkdirRobust([ savePath]);
    mkdirRobust([ savePath filesep 'XY'  ]);
    mkdirRobust([ savePath filesep 'ZY' ]);
    mkdirRobust([ savePath filesep 'ZX' ]);
    mkdirRobust([ savePath filesep 'three' ]);
end

XYFilesPattern=[savePath filesep 'XY' filesep 'XY_frame_nb%04d.png'];
YZFilesPattern=[savePath filesep 'ZY' filesep 'ZY_frame_nb%04d.png'];
XZFilesPattern=[savePath filesep 'ZX' filesep 'ZX_frame_nb%04d.png'];
ThreeFilesPattern=[savePath filesep 'three' filesep 'Three_frame_nb%04d.png'];

maxXBorder=MD.getDimensions('X');
maxYBorder=MD.getDimensions('Y');
maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
minXBorder=1;
minYBorder=1;
minZBorder=1;
frameNb=MD.nFrames_;
save([savePath filesep 'limits.mat'],'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');


process=ExternalProcess(MD,'printMIP');
process.setOutFilePaths({XYFilesPattern,YZFilesPattern,XZFilesPattern,[savePath filesep 'limits.mat'],ThreeFilesPattern});
MD.addProcess(process);
MD.save();


ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
minIntensityNorm=[];
maxIntensityNorm=[];
for chIdx=1:length(MD.channels_)
    vol=MD.getChannel(chIdx).loadStack(1);
    minIntensityNorm=[ minIntensityNorm min(vol(:))];
    maxIntensityNorm=[ maxIntensityNorm max(vol(:))];
end

parfor frameIdx=1:MD.nFrames_
%     progressText(frameIdx/MD.nFrames_,'Print MIP');
    maxXY=[];maxZY=[];maxZX=[];three=[];
    for chIdx=1:length(MD.channels_)
        vol=MD.getChannel(chIdx).loadStack(frameIdx);  
        [cmaxXY,cmaxZY,cmaxZX,cthree]=computeMIPs(vol,ZXRatio,minIntensityNorm(chIdx),maxIntensityNorm(chIdx));
        maxXY=[ maxXY cmaxXY];
        maxZY=[ maxZY cmaxZY];
        maxZX=[ maxZX cmaxZX];
        three=[ three cthree];        
    end
    
    % save the maximum intensity projections
    imwrite(maxXY, sprintfPath(XYFilesPattern,frameIdx), 'Compression', 'none');
    imwrite(maxZY, sprintfPath(YZFilesPattern,frameIdx), 'Compression', 'none');
    imwrite(maxZX, sprintfPath(XZFilesPattern,frameIdx), 'Compression', 'none');
    imwrite(three, sprintfPath(ThreeFilesPattern,frameIdx), 'Compression', 'none');
end

video = VideoWriter([savePath filesep 'three.avi']);
video.FrameRate = 4;  % Default 30
video.Quality = 90;    % Default 75

open(video)
for frameIdx=1:MD.nFrames_
    % save the maximum intensity projections
    three=imread(sprintfPath(ThreeFilesPattern,frameIdx));
    writeVideo(video,three)
%     fprintf('\b|\n');
end
close(video)