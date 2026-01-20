function out = pickReferenceLayerAndSave(md2D, beadChan, refMD, refBeadChan, varargin)
%pickReferenceLayerAndSave
% Given a 2D MovieData (single z-plane movie) and a reference z-stack MovieData,
% find the reference z-layer most similar to the bead image in md2D, save that
% layer as a reference image under md2D.outputDirectory_, and (optionally) try
% to register it to TFMPackage funParams (best-effort).
%
% out = pickReferenceLayerAndSave(mdZ1, beadChan, refMD, refBeadChan, ...)
% mdZ1: splitZMovie? ?? 2D MovieData (z=1)
% refMD: reference bead z-stack MovieData (zSize_ > 1)
% beadChan / refBeadChan: bead ?? index
% 
% out = pickReferenceLayerAndSave(mdZ1, beadChan, refMD, refBeadChan, ...
%     'targetFrame', 1, ...
%     'refFrame', 1, ...
%     'saveName', 'refBeads_fromStack.tif', ...
%     'registerTFM', true, ...
%     'verbose', true);
% 
% disp(out.bestZ)
% disp(out.savedPath)
%
% Name-Value options:
%   'targetFrame'   : frame index in md2D to use as template (default: 1)
%   'refFrame'      : frame index in refMD to scan across z (default: 1)
%   'cropROI'       : [x y w h] ROI applied to BOTH images before scoring (default: [])
%   'downsample'    : scalar >=1; if >1, images are resized by 1/downsample for scoring (default: 1)
%   'saveName'      : output filename (default: 'reference_beads.tif')
%   'saveFolder'    : subfolder under md2D.outputDirectory_ (default: 'reference')
%   'registerTFM'   : true/false, best-effort update TFMPackage funParams (default: true)
%   'verbose'       : true/false (default: true)
%
% Output struct out:
%   out.bestZ
%   out.scores
%   out.savedPath
%   out.templateFrame
%   out.refFrame

ip = inputParser;
ip.addRequired('md2D', @(x) isa(x,'MovieData'));
ip.addRequired('beadChan', @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addRequired('refMD', @(x) isa(x,'MovieData'));
ip.addRequired('refBeadChan', @(x) isnumeric(x) && isscalar(x) && x>=1);

ip.addParameter('targetFrame', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('refFrame', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('cropROI', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==4));
ip.addParameter('downsample', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('saveName', 'reference_beads.tif', @(x) ischar(x) || isstring(x));
ip.addParameter('saveFolder', 'reference', @(x) ischar(x) || isstring(x));
ip.addParameter('registerTFM', true, @(x) islogical(x) && isscalar(x));
ip.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
ip.parse(md2D, beadChan, refMD, refBeadChan, varargin{:});

tFrame   = ip.Results.targetFrame;
rFrame   = ip.Results.refFrame;
cropROI  = ip.Results.cropROI;
ds       = ip.Results.downsample;
saveName = char(ip.Results.saveName);
saveFolder = char(ip.Results.saveFolder);
doRegister = ip.Results.registerTFM;
verbose  = ip.Results.verbose;

% Basic checks
if isempty(refMD.zSize_) || refMD.zSize_ < 1
    error('refMD does not appear to be a z-stack (refMD.zSize_ empty or <1).');
end
if tFrame > md2D.nFrames_
    error('targetFrame exceeds md2D.nFrames_.');
end
if rFrame > refMD.nFrames_
    error('refFrame exceeds refMD.nFrames_.');
end

% Load template image from md2D
I0 = md2D.getChannel(beadChan).loadImage(tFrame);
I0 = toDoubleNorm(I0);

% Apply ROI if requested
if ~isempty(cropROI)
    I0 = imcrop(I0, cropROI);
end

% Downsample for scoring
if ds > 1
    I0s = imresize(I0, 1/ds);
else
    I0s = I0;
end

Z = refMD.zSize_;
scores = nan(Z,1);

if verbose
    fprintf('[pickReferenceLayerAndSave] Scanning %d z-layers...\n', Z);
end

for z = 1:Z
    J = refMD.getChannel(refBeadChan).loadImage(rFrame, z);
    J = toDoubleNorm(J);

    if ~isempty(cropROI)
        J = imcrop(J, cropROI);
    end

    % Match size (robust: crop to common min dims)
    [Icmp, Jcmp] = matchSize(I0s, (ds>1)*imresize(J,1/ds) + (ds==1)*J);

    % Similarity score (Pearson-like via corr2)
    % corr2 expects finite values; handle constant images
    if std(Icmp(:)) < eps || std(Jcmp(:)) < eps
        scores(z) = -Inf;
    else
        curScore2D = normxcorr2(Icmp, Jcmp);
        scores(z) = max(curScore2D(:));
    end
end

[~, bestZ] = max(scores);

% Load the best Z layer at full resolution (no ds), apply ROI, and save
bestImg = refMD.getChannel(refBeadChan).loadImage(rFrame, bestZ);

if ~isempty(cropROI)
    bestImg = imcrop(bestImg, cropROI);
end

outDir = fullfile(md2D.outputDirectory_, saveFolder);
if ~exist(outDir, 'dir'); mkdir(outDir); end
savedPath = fullfile(outDir, saveName);

% Preserve datatype as uint16 if possible
imwrite(castToUInt(bestImg), savedPath);

if verbose
    fprintf('[pickReferenceLayerAndSave] bestZ = %d, score = %.4f\n', bestZ, scores(bestZ));
    fprintf('[pickReferenceLayerAndSave] saved reference image: %s\n', savedPath);
end

% Best-effort: register to TFM package parameters if present
registerInfo = struct('attempted', false, 'updated', false, 'message', '');

if doRegister
    registerInfo.attempted = true;
    registerInfo = tryRegisterToTFM(md2D, savedPath, bestZ, registerInfo, verbose);
end

out = struct();
out.bestZ = bestZ;
out.scores = scores;
out.savedPath = savedPath;
out.templateFrame = tFrame;
out.refFrame = rFrame;
out.registerInfo = registerInfo;

end

%% -------- helper functions --------

function A = toDoubleNorm(A)
A = double(A);
% normalize roughly (avoid scale dominating corr if different gain)
A = A - mean(A(:));
sd = std(A(:));
if sd > 0, A = A / sd; end
end

function [A2, B2] = matchSize(A, B)
h = min(size(A,1), size(B,1));
w = min(size(A,2), size(B,2));
A2 = A(1:h,1:w);
B2 = B(1:h,1:w);
end

function U = castToUInt(I)
% If already integer, keep type; otherwise scale to uint16 safely
if isinteger(I)
    U = I;
else
    I = double(I);
    I = I - min(I(:));
    mx = max(I(:));
    if mx > 0, I = I / mx; end
    U = uint16(round(I * 65535));
end
end

function reg = tryRegisterToTFM(md2D, savedPath, bestZ, reg, verbose)
try
    if isempty(md2D.packages_)
        reg.message = 'md2D.packages_ is empty. Nothing to register.';
        return;
    end

    % Find a package whose class name contains 'TFM'
    pkgIdx = [];
    for k = 1:numel(md2D.packages_)
        if isempty(md2D.packages_{k}), continue; end
        cname = class(md2D.packages_{k});
        if contains(lower(cname), 'tfm')
            pkgIdx = k; break;
        end
    end
    if isempty(pkgIdx)
        reg.message = 'No TFMPackage found in md2D.packages_. Saved image only.';
        return;
    end

    pkg = md2D.packages_{pkgIdx};
    updated = false;

    % Best-effort: walk processes and patch funParams_ if fields exist
    for p = 1:numel(pkg.processes_)
        pr = pkg.processes_{p};
        if isempty(pr), continue; end

        if isprop(pr,'funParams_') && ~isempty(pr.funParams_) && isstruct(pr.funParams_)
            fp = pr.funParams_;

            % Common-ish field guesses:
            candidateFields = fieldnames(fp);
            % Set any obvious path-like fields
            for f = 1:numel(candidateFields)
                fn = candidateFields{f};
                low = lower(fn);
                if contains(low,'reference') && (contains(low,'path') || contains(low,'file'))
                    fp.(fn) = savedPath; updated = true;
                end
            end
            % Also try common names explicitly
            if isfield(fp,'referenceFramePath'), fp.referenceFramePath = savedPath; updated = true; end
            if isfield(fp,'refFramePath'),      fp.refFramePath      = savedPath; updated = true; end
            if isfield(fp,'referencePath'),     fp.referencePath     = savedPath; updated = true; end
            if isfield(fp,'referenceZ'),        fp.referenceZ        = bestZ;    updated = true; end

            if updated
                pr.funParams_ = fp; %#ok<NASGU>
                pkg.processes_{p} = pr;
            end
        end
    end

    if updated
        md2D.packages_{pkgIdx} = pkg;
        md2D.save();
        reg.updated = true;
        reg.message = 'Updated TFMPackage funParams_ (best-effort) and saved MovieData.';
        if verbose
            fprintf('[pickReferenceLayerAndSave] TFMPackage funParams updated (best-effort) and md2D saved.\n');
        end
    else
        reg.message = 'TFMPackage found but no matching funParams fields to update. Saved image only.';
    end

catch ME
    reg.message = sprintf('TFM registration attempt failed: %s', ME.message);
end
end