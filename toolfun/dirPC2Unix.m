function outImgDir = dirPC2Unix(inImgDir,imgDrive)
%dirPC2Unix: Convert image directories from PC format to Unix format.
%
% INPUT:
%    inImgDir: Image directories in PC format. It can be a cell array.
%    imgDrive: The disk drive name for the image directory in Unix format.
%
% See also: dirUnix2PC.

if ~strcmp(imgDrive(1),'/')
   imgDrive = ['/' imgDrive];
end

if ~iscell(inImgDir)
   outImgDir = {inImgDir};
else
   outImgDir = inImgDir;
end

for k = 1:length(outImgDir)
   fileSepInd = findstr('\',outImgDir{k});
   outImgDir{k}(fileSepInd) = '/';

   outImgDir{k} = [imgDrive outImgDir{k}(fileSepInd(1):end)];
end

if ~iscell(inImgDir)
   outImgDir = outImgDir{1};
end
