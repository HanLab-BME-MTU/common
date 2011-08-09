function stk2tiffDirs(varargin)
%
% function stk2tiffDirs(stkpath)
%
% Splits STKs in input directory into folders with TIFF files.
%
% Input: path to directory containing STK files
%
% Francois Aguet, 09/01/2010


ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('Crop', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.parse(varargin{:});
stkpath = ip.Results.path;
crop = ip.Results.Crop;

if isempty(stkpath)
   stkpath = uigetdir('Select directory containing the STK files:'); 
   if (stkpath == 0)
       return;
   end
end

% Recursive call if input is cell array
if iscell(stkpath), 
    cellfun(@(x)stk2tiffDirs(x,'Crop',crop),stkpath);
    return
end

stkpath = [stkpath filesep];
stkList = [dir([stkpath '*.tif']) dir([stkpath '*.tiff']) dir([stkpath '*.stk'])];

N = length(stkList);

for k = 1:N
    fprintf('Converting: %s\n', stkList(k).name);
    [~,stkname] = fileparts(stkList(k).name);
    dirname = [stkpath stkname filesep];
    [~,~] = mkdir(dirname);
    
    stack = stackRead([stkpath stkList(k).name]);
    
    if strcmpi(ip.Results.Crop, 'on')
       h = figure;
       imagesc(stack(:,:,1)); colormap(gray(256)); axis image;
       hr = imrect();
       pos = round(wait(hr));
       close(h);
       stack = stack(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
    end
    
    nf = size(stack,3);  
    for z = 1:nf
        imwrite(stack(:,:,z), [dirname stkname '_' num2str(z, ['%0' num2str(length(num2str(nf))) 'd']) '.tif'], 'tif');
    end;
end;