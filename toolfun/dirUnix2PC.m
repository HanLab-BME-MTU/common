function outImgDir = dirUnix2PC(inImgDir,imgDrive)
%dirPC2Unix: Convert image directories from Unix format to PC format.
%
% INPUT:
%    inImgDir: Image directories in PC format. It can be a cell array.
%    imgDrive: The disk drive letter for the image directory in PC format.
% 
% See also: dirPC2Unix.

colonInd = findstr(':',imgDrive);
if isempty(colonInd)
   imgDrive = [imgDrive ':'];
end

if ~iscell(inImgDir)
   outImgDir = {inImgDir};
else
   outImgDir = inImgDir;
end

for k = 1:length(outImgDir)
   fileSepInd = findstr('/',outImgDir{k});
   outImgDir{k}(fileSepInd) = '\';

   mntInd = findstr('/mnt/',outImgDir{k});
   if isempty(mntInd)
      outImgDir{k} = [imgDrive outImgDir{k}(fileSepInd(2):end)];
   else
      outImgDir{k} = [imgDrive outImgDir{k}(fileSepInd(3):end)];
   end
end

if ~iscell(inImgDir)
   outImgDir = outImgDir{1};
end
