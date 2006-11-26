function mkRGDuoChannelImg(varargin)
%mkRGDuoChannelImg: This function creates Red/Green duo-channel images.
%
% SYNOPSIS:
%    mkRGDuoChannelImg;
%    mkRGDuoChannelImg('outFileName',outFileName,'redImgIScale',redImgIScale,'greenImgIScale',greenImgIScale);
%
% OPTIONAL PAR/VALUE INPUT:
%        PAR
%    -------------
%    'outFileName': The common prefix name for the output files. Default is the
%                   same name as the image file name of the first channel with 'rg_' as prefix.
%    'redImgIScale': (red channel) A vector of two numbers between 0 and 1 that specifies
%                    the percentage of cut-offs at the lower and upper end of image intensity.
%    'greenImgIScale': (green channel) A vector of two numbers between 0 and 1 that specifies
%                    the percentage of cut-offs at the lower and upper end of image intensity.
%
% AUTHOR: Lin Ji, Nov. 26, 2005

[redFirstImgFile redInDir filterIndex] = uigetfile('*.*','Pick the first image for red channel');
if isequal(redFirstImgFile,0) || isequal(redInDir,0)
   disp('User pressed cancel. No file is selected.');
   return;
end

[path,redFirstImgFName,no,ext] = getFilenameBody(redFirstImgFile);
redFirstImgFIndex = str2num(no);

%Default output file name.
outFileName = ['rg_' redFirstImgFName];

[redImgFList,index] = getNamedFiles(redInDir,redFirstImgFName);
selIndex = find(index>=redFirstImgFIndex);
redImgFIndex    = index(selIndex);
redImgFList = redImgFList(selIndex);

numRedImgFiles  = length(redImgFList);

[greenFirstImgFile greenInDir filterIndex] = uigetfile('*.*','Pick the first image for green channel');
if isequal(greenFirstImgFile,0) || isequal(greenInDir,0)
   disp('User pressed cancel. No file is selected.');
   return;
end

[path,greenFirstImgFName,no,ext] = getFilenameBody(greenFirstImgFile);
greenFirstImgFIndex = str2num(no);

[greenImgFList,index] = getNamedFiles(greenInDir,greenFirstImgFName);
selIndex = find(index>=greenFirstImgFIndex);
greenImgFIndex    = index(selIndex);
greenImgFList = greenImgFList(selIndex);

numGreenImgFiles  = length(greenImgFList);

if numRedImgFiles ~= numGreenImgFiles
   error('The number of image files in red and green channels do not match.');
end

parentDir = [redInDir filesep '..'];
outDir = uigetdir(parentDir,'Select output directory');

if isequal(outDir,0)
   disp('User pressed cancel. No output directory is selected.');
   return;
end

%Default
redImgIScale   = [0 1];
greenImgIScale = [0 1];

if nargin > 0
   for kk = 1:2:nargin
      switch varargin{kk}
         case 'outFileName'
            outFileName = varargin{kk+1};
         case 'redImgIScale'
            redImgIScale = varargin{kk+1};
         case 'greenImgIScale'
            greenImgIScale = varargin{kk+1};
      end
   end
end

figH = figure; hold off;
for kk = 1:numRedImgFiles
   [path,imgFName,no,ext] = getFilenameBody(redImgFList{kk});
   imgIndexStr = no;

   redImg = double(imread([redInDir filesep redImgFList{kk}]));
   redImg = redImg-min(redImg(:));
   rgImg = zeros([size(redImg) 3]);

   maxRedImgI = max(redImg(:));

   redImgF = redImg;
   redImgF(find(redImg>redImgIScale(2)*maxRedImgI)) = redImgIScale(2)*maxRedImgI;
   redImgF(find(redImg<redImgIScale(1)*maxRedImgI)) = redImgIScale(1)*maxRedImgI;

   maxRedImgI = max(redImgF(:));
   minRedImgI = min(redImgF(:));

   redImgF = (redImgF-minRedImgI)/(maxRedImgI-minRedImgI);

   greenImg = double(imread([greenInDir filesep greenImgFList{kk}]));
   greenImg = greenImg-min(greenImg(:));
   rgImg = zeros([size(greenImg) 3]);

   maxGreenImgI = max(greenImg(:));

   greenImgF = greenImg;
   greenImgF(find(greenImg>greenImgIScale(2)*maxGreenImgI)) = greenImgIScale(2)*maxGreenImgI;
   greenImgF(find(greenImg<greenImgIScale(1)*maxGreenImgI)) = greenImgIScale(1)*maxGreenImgI;

   maxGreenImgI = max(greenImgF(:));
   minGreenImgI = min(greenImgF(:));

   greenImgF = (greenImgF-minGreenImgI)/(maxGreenImgI-minGreenImgI);

   rgImg(:,:,1) = redImgF;
   rgImg(:,:,2) = greenImgF;

   figure(figH); hold off;
   imshow(rgImg,[]);

   rgImgFile = [outDir filesep outFileName imgIndexStr '.fig'];
   saveas(figH,rgImgFile,'fig');
end
