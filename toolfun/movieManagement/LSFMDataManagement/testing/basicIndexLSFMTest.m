root='/project/bioinformatics/Danuser_lab/shared/proudot/dataManagement/exampleDataSets/Actin_BF_noBF_mix/'

% Single movie

%% test useBF
file=[root 'Cell2/1_CAM01_000000.tif'];
ML=indexLSFMData(file,fileparts(file),'useBF',true)
MD=ML.getMovie(1);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.nFrames_==GT.MD.nFrames_)
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)

%% test useBF channel and regexp no MIP
fileReg=[root 'Cell2/1_CAM0{ch}_*.tif']
ML=indexLSFMData(fileReg,fileparts(fileReg),'useBF',true,'createMIP',false)
MD=ML.getMovie(1);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.nFrames_==GT.MD.nFrames_)
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)


%% test channel based building
fileReg=[root 'Cell2/1_CAM0{ch}_*.tif']
ML=indexLSFMData(fileReg,fileparts(fileReg),'useBF',false,'copyFile',true,'createMIP',false,'lateralPixelSize',GT.MD.pixelSize_,'axialPixelSize',GT.MD.pixelSizeZ_)
MD=ML.getMovie(1);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.nFrames_==GT.MD.nFrames_)
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)

% Multiple movie

%% test channel based building on multiple movies
fileReg=[root 'Cell*/1_CAM0{ch}_*']
ML=indexLSFMData(fileReg,root,'useBF',false,'copyFile',true,'createMIP',false,'lateralPixelSize',GT.MD.pixelSize_,'axialPixelSize',GT.MD.pixelSizeZ_)
MD=ML.getMovie(2);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)

%% test BF based building on multiple movies
fileReg=[root 'Cell*/1_CAM0{ch}_*']
ML=indexLSFMData(fileReg,root,'useBF',true,'copyFile',true,'createMIP',false)
MD=ML.getMovie(2);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)

%% test BF based building on multiple movies no root
fileReg=[root 'Cell*/1_CAM01_000000.tif']
ML=indexLSFMData(fileReg,'useBF',true,'copyFile',true,'createMIP',false)
MD=ML.getMovie(2);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)

%% test BF based building on multiple movies no root
fileReg=[root 'Cell*/1_CAM0{ch}_*']
ML=indexLSFMData(fileReg,'useBF',true)
MD=ML.getMovie(2);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.pixelSize_==GT.MD.pixelSize_)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==GT.MD.pixelSizeZ_)


break; 

%% WARNING: file moved after that test
%% test channel based building just optional argument
fileReg=[root 'Cell2/1_CAM0{ch}_*.tif']
ML=indexLSFMData(fileReg,fileparts(fileReg))
MD=ML.getMovie(1);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.nFrames_==GT.MD.nFrames_)
assert(MD.pixelSize_==1)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==1)

%% test channel based building just required argument
fileReg=[root 'Cell2/1_CAM0{ch}_*.tif']
ML=indexLSFMData(fileReg)
MD=ML.getMovie(1);
GT=load([root 'analysis/cell4-MovieData-GT.mat']);
assert(MD.nFrames_==GT.MD.nFrames_)
assert(MD.pixelSize_==1)
assert(MD.zSize_==GT.MD.zSize_)
assert(MD.pixelSizeZ_==1)
