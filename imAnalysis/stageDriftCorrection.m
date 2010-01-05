function T = stageDriftCorrection(inputFileList, sigmaPSF)
% stageDriftCorrection returns an array containing all drifts between each
% consecutive pair of images.
%
% SYNOPSIS T = stageDriftCorrection(...)
%
% INPUT    inputFileList: cell array of all image filenames (including full
%                         path name).
%
%          sigmaPSF: hald-with of the point spread function (standard
%          deviation of a Gaussian model PSF).
%
% OUTPUT   T: an array of the same size of inputFileList minus 1 containing
%          all 2D drifts between images.
%
% Sylvain Berlemont, December 12th, 2009

T = [];

if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([filename filesep pathname]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if numel(inputFileList) < 2
    error('Number of images is less than 2.');
end

n = numel(inputFileList);
pts = cell(n, 1);
T = zeros(n - 1, 2);

%
% Step 1: Subpixel spot detection
%

for i = 1:n
    % Read image.
    I = double(imread(inputFileList{i}));
    
    % Denoise image using Wavelet A Trou.
    Irec = awtDenoising(I, [], 0, 3);
    % Find local maximum
    w = 2 * ceil(sigmaPSF) + 1;
    locMax = find(locmax2d(Irec, [w, w]));

    if isempty(locMax)
        error('Frame %d does not contain any point.', ind);
    end
    
    [y x] = ind2sub(size(Irec), locMax);
    
    % Subpixel detection
    estimates = fitMixModel(I, [y x], sigmaPSF, I(locMax), ...
        mean(nonzeros(I - Irec)));
    
    pts{i} = estimates(:,1:2);
end

%
% Step 2: Point registration
%

numIter = 10;
tol = 1e-4;
tolR = 1e-4;

for i = 1:n-1
    n1 = size(pts{i}, 1);
    X1 = [pts{i}(:, 1:2) zeros(n1, 1)];
    n2 = size(pts{i+1}, 1);
    X2 = [pts{i+1}(:, 1:2) zeros(n2, 1)];
    
    if n1 > n2
        [Ri, Ti] = computeICP(X1, X2, numIter, tol);
        Ti = -Ti;
    else
        [Ri, Ti] = computeICP(X2, X1, numIter, tol);
    end
    
    Ri = round(Ri / tolR) * tolR;
    
    if ~all(all(Ri == eye(3)))
        warning('Warning: significant rotation effect detected between frame %d-%d.',...
            i, i+1); %#ok<WNTAG>
        T(i, :) = NaN;
    else
        T(i, :) = Ti(1:2);
    end
    % Free allocated memory for the kd-tree.
    kdtree([], [], treeRoot);
end