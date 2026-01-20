function mdList = splitZMovie(movieData, outRoot, varargin)
%SPLITZMOVIE Split a z-stack MovieData into separate 2D MovieData objects (one per z-plane).
%
% mdList = splitZMovie(movieData, outRoot)
% mdList = splitZMovie(movieData, outRoot, 'zList', 1:3, 'cropROI', [x y w h], 'cropTOI', 1:T)
%
% Inputs
%   movieData : MovieData (must have zSize_ >= 1)
%   outRoot   : root output directory. Each z-plane will be saved under:
%               outRoot/Z01, outRoot/Z02, ...
%
% Name-Value pairs
%   'zList'   : vector of z indices to export (default: 1:movieData.zSize_)
%   'cropROI' : [xmin ymin width height] (default: full image)
%   'cropTOI' : vector of time indices (default: 1:movieData.nFrames_)
%   'prefix'  : folder prefix (default: 'Z')
%   'digits'  : number of digits in folder numbering (default: 2)
%   'verbose' : true/false (default: true)
%
% Output
%   mdList    : cell array of MovieData objects, one per z-plane

ip = inputParser;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addRequired('outRoot', @(x) ischar(x) || isstring(x));

defaultROI = [1 1 movieData.imSize_(end:-1:1)];
defaultTOI = 1:movieData.nFrames_;

ip.addParameter('zList', [], @(x) isnumeric(x) && isvector(x));
ip.addParameter('cropROI', defaultROI, @(x) isnumeric(x) && isvector(x) && numel(x)==4);
ip.addParameter('cropTOI', defaultTOI, @(x) isnumeric(x) && isvector(x));
ip.addParameter('prefix', 'Z', @(x) ischar(x) || isstring(x));
ip.addParameter('digits', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
ip.parse(movieData, outRoot, varargin{:});

outRoot = char(ip.Results.outRoot);
cropROI = ip.Results.cropROI;
cropTOI = ip.Results.cropTOI;
prefix = char(ip.Results.prefix);
digits = ip.Results.digits;
verbose = ip.Results.verbose;

% Validate zSize
if isempty(movieData.zSize_) || movieData.zSize_ < 1
    error('splitZMovie:NoZ', 'movieData.zSize_ is empty or < 1. This MovieData does not appear to be a z-stack.');
end

% Default zList
zList = ip.Results.zList;
if isempty(zList)
    zList = 1:movieData.zSize_;
end
zList = unique(zList(:)'); % row
if any(zList < 1) || any(zList ~= round(zList))
    error('splitZMovie:BadZList', 'zList must contain positive integer z-indices.');
end
if max(zList) > movieData.zSize_
    error('splitZMovie:ZOutOfRange', 'Requested zList exceeds movieData.zSize_ (%d).', movieData.zSize_);
end

% Create root dir
if ~exist(outRoot, 'dir')
    mkdir(outRoot);
end

mdList = cell(1, numel(zList));

for k = 1:numel(zList)
    z = zList(k);
    zTag = sprintf(['%s%0' num2str(digits) 'd'], prefix, z);
    zOut = fullfile(outRoot, zTag);
    mkClrDir(zOut)

    if verbose
        fprintf('[splitZMovie] Exporting z=%d -> %s\n', z, zOut);
    end

    % Call cropMovie with cropZOI
    mdList{k} = cropMovie(movieData, zOut, ...
        'cropROI', cropROI, ...
        'cropTOI', cropTOI, ...
        'cropZOI', z);
end

end