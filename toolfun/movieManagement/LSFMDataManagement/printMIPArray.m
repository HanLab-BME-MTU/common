function printMIPArray(ML,varargin)
% Print an array of XY MIP for all Movie in the MovieList. The resulting
% images and movies are saved in the movieList output directory.

ip = inputParser;
ip.addRequired('ML',@(MD) isa(ML,'MovieList'));
ip.addParamValue('maxWidth',1600,@isnumeric);
ip.addParamValue('maxHeight',900,@isnumeric);
ip.addParamValue('MIPSize',250,@isnumeric);
ip.addParamValue('invertLUT',0,@islogical);
ip.parse(ML,varargin{:});
p=ip.Results;
% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');


% set parameters
stripeSize = 8; % the width of the stripes in the image that combines all three maximum intensity projections
%the stripe color, a number between 0 (black) and 1 (white).  (If you're not using all of the bit depth then white might be much lower than 1, and black might be higher than 0.)
% if(p.invertLUT); 
%     stripeColor = 1; 
% else
    stripeColor = 0; 
% end;

savePath=[ML.outputDirectory_ filesep 'MIPArray'];
if ~isdir(savePath) 
    mkdir(savePath)
end

MD=ML.getMovie(1);
nameCells=MD.getChannel(1).getImageFileNames;
% fprintf('printing MIP:');
% fprintf(['\n' repmat('.',1,MD.nFrames_) '\n\n']);
minIntensityNorm=[];
maxIntensityNorm=[];
for chIdx=1:length(MD.channels_)
    vol=MD.getChannel(chIdx).loadStack(1);
    minIntensityNorm=[ minIntensityNorm min(vol(:))];
    maxIntensityNorm=[ maxIntensityNorm max(vol(:))];
end
        

video = VideoWriter([savePath filesep 'MIPArray.avi']);
myVideo.FrameRate = 4;  % Default 30
myVideo.Quality = 90;    % Default 75
open(video)

% Collect  info to build array and intensity normalization
% MDValid=ones(1,ML.getSize());
frameNb=nan(1,ML.getSize());
chNb=nan(1,ML.getSize());
chMinIntensityNorm=cell(1,ML.getSize());
chMaxIntensityNorm=cell(1,ML.getSize());
for MDIdx=1:ML.getSize
    MD=ML.getMovie(MDIdx);
%     if(MD.zSize_==1)
%         warning('Movie ' num2str(MDIdx) ' is a 2D movie, No MIP produced.');
%         MDValid(MDIdx)=0;
%     end
    frameNb(MDIdx)=MD.nFrames_;
    chNb(MDIdx)=length(MD.channels_);
    chMinIntensityNorm{MDIdx}=[];
    chMaxIntensityNorm{MDIdx}=[];
    for chIdx=1:length(MD.channels_)
        vol=MD.getChannel(chIdx).loadStack(1);
        chMinIntensityNorm{MDIdx}=[ chMinIntensityNorm{MDIdx} min(vol(:))];
        chMaxIntensityNorm{MDIdx}=[ chMaxIntensityNorm{MDIdx} max(vol(:))];
    end
end


%% Build array for movie location
maxMoviePerLine=floor(p.maxWidth/(p.MIPSize*max(chNb)+stripeSize));
maxColumns=floor(p.maxHeight/(p.MIPSize+stripeSize));
MDArray=(1:ML.getSize);
untruncatedArraySize=[ceil(ML.getSize()/maxMoviePerLine), min(ML.getSize(),maxMoviePerLine)];
if(untruncatedArraySize(1)*untruncatedArraySize(2)>ML.getSize)
    MDArray(untruncatedArraySize(1)*untruncatedArraySize(2))=0;
end
MDArray=reshape(MDArray,untruncatedArraySize);
MDArray=MDArray(1:min(end,maxColumns),:);

%%
for frameIdx=1:max(frameNb)
    maxXY=[];
    for i=1:size(MDArray,1)
        movieLine=[];
        for MDIdx=MDArray(i,MDArray(i,:)>0)
            MD=ML.getMovie(MDIdx);
            mMaxXY=[];
            ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
            for chIdx=1:length(MD.channels_)
                vol=MD.getChannel(chIdx).loadStack(min(frameIdx,MD.nFrames_));
                [cmaxXY]=computeMIPs(vol,ZXRatio,chMinIntensityNorm{MDIdx}(chIdx),chMaxIntensityNorm{MDIdx}(chIdx));
                cmaxXY=imresize(cmaxXY,p.MIPSize/max(size(cmaxXY)));
                cmaxXY=padarray(cmaxXY,[max(0,p.MIPSize-size(cmaxXY,1)) 0],'post');
                mMaxXY=[ mMaxXY cmaxXY];            
            end
            movieLine=[movieLine mMaxXY];
        end
        %%
        movieLine=padarray(movieLine,[0 p.maxWidth-size(movieLine,2)],'post');
        %%
        maxXY=[maxXY; stripeColor*ones(stripeSize,p.maxWidth); movieLine];
    end

    if(p.invertLUT)
        maxXY = imcomplement(maxXY);
    end
    
    imwrite(maxXY, [savePath filesep 'MIPArray_' num2str(frameIdx,'%04d') '.tif' ], 'Compression', 'none');
    writeVideo(video,maxXY)
    
end
close(video)

function mkdir(path)
system(['mkdir ' path]);


