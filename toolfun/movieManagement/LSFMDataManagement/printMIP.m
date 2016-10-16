function printMIP(MD)
% Title says it all -- PR, augmented from Meghan D. 2015

ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.parse(MD);



% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

if(MD.zSize_==1)
    error('This seems to be a 2D movie, No MIP produced.');
end

savePath=[MD.outputDirectory_ filesep 'MIP'];
if ~isdir(savePath) || ~isdir([savePath filesep 'XY']) || ~isdir([savePath filesep 'ZY']) || ~isdir([savePath filesep 'ZX']) || ~isdir([savePath filesep 'Three'])
    mkdir(savePath)
    mkdir([savePath filesep 'XY'])
    mkdir([savePath filesep 'ZY'])
    mkdir([savePath filesep 'ZX'])
    mkdir([savePath filesep 'Three'])
end
nameCells=MD.getChannel(1).getImageFileNames;
% fprintf('printing MIP:');
% fprintf(['\n' repmat('.',1,MD.nFrames_) '\n\n']);
ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
minIntensityNorm=[];
maxIntensityNorm=[];
for chIdx=1:length(MD.channels_)
    vol=MD.getChannel(chIdx).loadStack(1);
    minIntensityNorm=[ minIntensityNorm min(vol(:))];
    maxIntensityNorm=[ maxIntensityNorm max(vol(:))];
end

threeVideo = VideoWriter([savePath filesep 'three.avi']);
myVideo.FrameRate = 4;  % Default 30
myVideo.Quality = 90;    % Default 75
open(threeVideo)

for frameIdx=1:MD.nFrames_
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
    imwrite(maxXY, [savePath filesep 'XY' filesep 'XY_' nameCells{frameIdx} ], 'Compression', 'none');
    imwrite(maxZY, [savePath filesep 'ZY' filesep 'ZY_'  nameCells{frameIdx}], 'Compression', 'none');
    imwrite(maxZX, [savePath filesep 'ZX' filesep 'ZX_'  nameCells{frameIdx}], 'Compression', 'none');
    imwrite(three, [savePath filesep 'Three' filesep nameCells{frameIdx}], 'Compression', 'none');
    writeVideo(threeVideo,three)
%     fprintf('\b|\n');
end
close(threeVideo)

