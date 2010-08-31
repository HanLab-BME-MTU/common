function stk2tiffDirs(stkpath)
%
% function stk2tiffDirs(stkpath)
%
% Splits STKs in input directory into folders with TIFF files.
%
% Input: path to directory containing STK files
%
% Francois Aguet, 01/20/2010

if (nargin == 0 || isempty(stkpath))
   stkpath = uigetdir('Select directory containing the STK files:'); 
   if (stkpath == 0)
       return;
   end;
end;

stkpath = [stkpath filesep];
stkList = [dir([stkpath '*.tif']) dir([stkpath '*.stk'])];

N = length(stkList);

for k = 1:N
    fprintf('Converting: %s\n', stkList(k).name);
    stkname = strrep(stkList(k).name(1:end-4), ' ', '_');
    dirname = [stkpath stkname filesep];
    if ~(exist(dirname, 'dir')==7)
        mkdir(dirname);
    end;
    stack = stackRead([stkpath stkList(k).name]);
    nf = size(stack,3);    
    for z = 1:nf
        imwrite(stack(:,:,z), [dirname stkname '_' num2str(z, ['%0' num2str(length(num2str(nf))) 'd']) '.tif'], 'tif');
    end;
end;