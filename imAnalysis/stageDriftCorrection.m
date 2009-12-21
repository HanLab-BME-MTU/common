function T = stageDriftCorrection(inputFileList, sigmaPSF, numIter, tol)
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
%          numIter: maximum number of iterations for the registration
%          process (Iterative Closest Point).
%
%          tol: maximum value |R - Id| underwhich any rotation R between the
%          sets of points is considered to be unsignificant. For a standard
%          microscope stage, only translation should be expected.
%
% OUTPUT   T: an array of the same size of inputFileList minus 1 containing
%          all 2D drifts between images.
%
% DEPENDENCES: extern/icp/*
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

% Number of iterations for ICP.
if nargin < 2 || ~isinteger(numIter)
    numIter = 50;
end

% Tolerance
if nargin < 3 || ~isnumeric(tol)
    tol = 1e-4;
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

for i = 1:n-1
    % Define model
    n_model = size(pts{i}, 1);
    model = [pts{i}(:, 2), pts{i}(:, 1) zeros(n_model, 1)]';
    % Define data
    n_data = size(pts{i+1}, 1);
    data = [pts{i+1}(:, 2), pts{i+1}(:, 1) zeros(n_data, 1)]';
    % weights
    weights = ones(1, n_data);
    % randvec
    rndvec = uint32(randperm(n_data)-1);
    % sizerand, number of point-matchings in each iteration.
    sizernd = ceil(1.45*n_data); % = min(n_model, n_data);
    % Create the kd-tree, TreeRoot is the pointer to the kd-tree
    [~, ~, treeRoot] = kdtree(model', []);
    % Run the ICP algorithm.
    [Ri, Ti] = icpCpp(model, data, weights, rndvec, sizernd, treeRoot, numIter);
    if norm(Ri(1:2, 1:2) - eye(2)) > tol
        warning('Warning: significant rotation effect detected between frame %d-%d.',...
            i, i+1); %#ok<WNTAG>
        T(i, :) = NaN;
    else
        T(i, :) = Ti(1:2);
    end
    % Free allocated memory for the kd-tree.
    kdtree([], [], treeRoot);
end