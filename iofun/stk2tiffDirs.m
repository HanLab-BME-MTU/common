function stk2tiffDirs(varargin)
% stktiffDirs splits STKs in input directory into folders with TIFF files.
% 
% Synopsis:    stk2tiffDirs(path)
%              stk2tiffDirs(path,'Crop','on')
%
% Input:
%      path : optional - path to the directory containing STK files. Can be
%      a string or a cell array of strings. If not input, the user will be
%      asked to select a folder.
%
%      Optional parameter/value pairs
%           Crop ('on'/'off') :  if 'on', a window will open asking the
%           user to select the region to crop.
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