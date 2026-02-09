function T = analyze_immunostain(opts)
% ANALYZE_IMMUNOSTAIN ? robust immunofluorescence quantification pipeline
% Quantifies Nucleus, Cytoplasm (nucleus-excluded), and Whole-Cell
% intensities for ALL channels in an OME-TIFF (multi-series = XY positions).
%
% SAVES (under <prefix>_Results):
%   Nucleus/<prefix>_Nucleus.csv
%   Cytoplasm/<prefix>_Cytoplasm.csv
%   WholeCell/<prefix>_WholeCell.csv
%   Masks/<prefix>_posXX_{nucleus|cytoplasm|wholecell}.png
%   <prefix>_ALL.csv (unified table with Compartment column)
%   Logs/method_log.txt (choices per position, pixel size, crop, background)
%   (optional) BgSub/<prefix>_posXX_<chan>.tif ? bg-subtracted channels
%   NEW:
%     - Masks/<prefix>_posXX_raw_<chan>.png (cropped raw channels)
%     - Masks/<prefix>_posXX_overlay.png (segmentation overlay)
%
% INTERACTION (no sliders):
%   Nucleus: shows 3 choices (Otsu / Otsu×0.8 / Adaptive(0.5)) ? locked to Otsu here
%   Cytoplasm: distance-based expansion seeded at nuclei, constrained by intensity
%  REQUIRED: pick a background image (single or multi-channel) ? no skip
%  Optional: global crop (integer indexing; no mismatch)
%  Optional: auto-detected channel names; can be renamed
%
% REQUIREMENTS:
%    bfmatlab (Bio-Formats): bfGetReader, bfGetPlane, etc.
%
% RETURN:
%   T : master table (Nucleus + Cytoplasm + Whole-Cell) with Compartment column.
%
% Author:

if nargin < 1, opts = struct; end
opts = fillDefault(opts, getDefaultOpts());

%% -------- File selection ----------
[fn, fp] = uigetfile({'*.ome.tif;*.ome.tiff;*.tif;*.tiff','OME-TIFF / TIF'}, 'Select OME-TIFF');
if isequal(fn,0), T = table(); disp('Canceled.'); return; end
imgPath = fullfile(fp, fn);

%% -------- Output prefix ----------
if isempty(opts.outputPrefix)
    [~, core, ~] = fileparts(fn);
    def = core;
else
    def = opts.outputPrefix;
end
if opts.askOutputPrefix
    answ = inputdlg({'Output prefix (used for folders & files):'}, 'Output Prefix', 1, {def});
    if ~isempty(answ) && ~isempty(answ{1}), opts.outputPrefix = answ{1}; else, opts.outputPrefix = def; end
end

%% -------- Output folders ----------
outRoot = fullfile(fp, [opts.outputPrefix '_Results']);
dirs = {'WholeCell','Nucleus','Cytoplasm','Masks','Logs'};
if opts.saveBgSubChannels, dirs{end+1} = 'BgSub'; end
for k = 1:numel(dirs)
    outDir.(dirs{k}) = fullfile(outRoot, dirs{k}); %#ok<AGROW>
    if ~exist(outDir.(dirs{k}),'dir'), mkdir(outDir.(dirs{k})); end
end

%% -------- Reader & pixel size ----------
reader = bfGetReader(imgPath);
nSeries = reader.getSeriesCount();

if isempty(opts.pixelSize_um)
    try
        px = reader.getMetadataStore().getPixelsPhysicalSizeX(0);
        if ~isempty(px), opts.pixelSize_um = double(px.value()); end
    catch
    end
end
if isempty(opts.pixelSize_um), opts.pixelSize_um = 0.108; end  % fallback if metadata missing
pxArea_um2 = opts.pixelSize_um^2;

%% -------- Channel names (auto + optional rename) ----------
nC = reader.getSizeC();
chanNames = cell(1,nC);
for c = 1:nC
    chanNames{c} = sprintf('Ch%d', c); % fallback
    try
        em = reader.getMetadataStore().getChannelEmissionWavelength(0,c-1);
        ex = [];
        if isempty(em), ex = reader.getMetadataStore().getChannelExcitationWavelength(0,c-1); end
        if ~isempty(em)
            chanNames{c} = sprintf('Ch%d', round(double(em.value())));
        elseif ~isempty(ex)
            chanNames{c} = sprintf('Ch%d', round(double(ex.value())));
        end
    catch
    end
end

if opts.askRenameChannels
    prompt = strcat("Enter name for ", chanNames, " (leave blank to keep)");
    answ = inputdlg(prompt, 'Channel Names', 1, chanNames);
    if ~isempty(answ)
        for c = 1:min(nC,numel(answ))
            if ~isempty(answ{c}), chanNames{c} = matlab.lang.makeValidName(answ{c}); end
        end
    end
end

%% -------- REQUIRED background selection (no skipping) ----------
[bgfn, bgfp] = uigetfile({'*.ome.tif;*.ome.tiff;*.tif;*.tiff','Background image'}, ...
                          'Select background image (single or multi-channel)');
if isequal(bgfn,0)
    try, reader.close(); catch, end
    error('Background image is required. Aborting.');
end
bgPath = fullfile(bgfp, bgfn);
bgMean = computeBackgroundMeans_fromImage(bgPath);

%% -------- Optional global crop (robust integer indexing) ----------
cropRect = [];
if opts.askCrop
    reader.setSeries(0);
    [stackC, ~] = readMultichannel(reader, 1);
    chShow = stackC(:,:,min(2, size(stackC,3)));
    f = figure('Name','Crop (Optional)');
    imshow(chShow,[]);
    title('Drag rectangle to crop; double-click to confirm; press Enter to skip');
    rect = getrect();
    if ishandle(f), close(f); end
    if ~isempty(rect) && rect(3)>1 && rect(4)>1
        cropRect = round(rect);
        fprintf('Using crop rect (rounded): [x=%d,y=%d,w=%d,h=%d]\n', cropRect);
    end
end

%% -------- Master tables & logs ----------
T_cell = table(); T_nuc = table(); T_cyto = table();
logLines = {};

%% -------- Iterate series ----------
for s = 1:nSeries
    reader.setSeries(s-1);
    [stack, ~] = readMultichannel(reader, 1);

    % Apply crop (if any)
    if ~isempty(cropRect)
        stack = cropStackByRect(stack, cropRect);
    end

    % ---- NEW: Save RAW (cropped) channels as PNG ----
    for c = 1:nC
        rawPath = fullfile(outDir.Masks, sprintf('%s_pos%02d_raw_%s.png', ...
            opts.outputPrefix, s, chanNames{c}));
        imwrite_uint16(rawPath, stack(:,:,c));
    end

    % ---- Background subtraction (image-based only, no rolling-ball) ----
    [stack, bgMethod, bgValue, bgReason] = enhancedBackgroundSubtraction(stack, bgMean, chanNames);

    % Save bg-subtracted channels if requested
    if opts.saveBgSubChannels
        for c = 1:nC
            imwrite_uint16(fullfile(outDir.BgSub, ...
                sprintf('%s_pos%02d_%s.tif', opts.outputPrefix, s, chanNames{c})), ...
                stack(:,:,c));
        end
    end

    % Split channels into cell array (for segmentation & intensity)
    channels = cell(1,nC);
    for c = 1:nC, channels{c} = stack(:,:,c); end

    % Nucleus channel
    nucIdx = 1;
    if isfield(opts.channelNames,'nuc405'), nucIdx = opts.channelNames.nuc405; end
    nucIdx = min(max(1, nucIdx), nC);
    chNuc = channels{nucIdx};

    % Cytoplasm reference channel
    if isfield(opts.channelNames,'piezo488')
        cytoRef = min(max(1,opts.channelNames.piezo488), nC);
    elseif nC >= 2
        cytoRef = 2;
    else
        cytoRef = 1;
    end

    satWarn(chNuc, chanNames{nucIdx});
    satWarn(channels{cytoRef}, chanNames{cytoRef});

    % Nucleus segmentation
    I = imgaussfilt(chNuc, opts.gaussSigma);
    [nucMask, nucLabels, thrChoice, thrInfo] = segmentNucleiInteractive(I, opts, s, pxArea_um2);

    if opts.useWatershed
        % Watershed used only for splitting touching nuclei
        nucLabels = splitAfterThreshold(nucMask, opts, pxArea_um2);
        nucMask   = nucLabels > 0;
    else
        [nucLabels, ~] = bwlabel(nucMask);
    end

    % Cytoplasm segmentation
    J = imgaussfilt(channels{cytoRef}, 1);
    [cytMask, cytLabels, sensUsed] = cytoUI_distanceConstrained(J, nucLabels, opts, s);

    % Whole-cell labels
    cellMask = nucMask | cytMask; %#ok<NASGU>
    cellLabels = nucLabels;
    m = cytLabels>0;
    cellLabels(m) = cytLabels(m);

    stats = regionprops(nucLabels,'PixelIdxList');
    nCells = numel(stats);

    colsMeta = {'Area_um2','Position','NucleusID','ThrChoice','ThrInfo','CytSens'};
    Tn = table('Size',[nCells numel(colsMeta)], ...
        'VariableTypes',{'double','double','double','string','string','double'}, ...
        'VariableNames',colsMeta);
    Tc = Tn; Tw = Tn;

    for i = 1:nCells
        nucPix  = stats(i).PixelIdxList;
        cytoPix = find(cytLabels == i);
        cellPix = find(cellLabels == i);

        % Areas
        Tn.Area_um2(i) = numel(nucPix)  * pxArea_um2;
        Tc.Area_um2(i) = numel(cytoPix) * pxArea_um2;
        Tw.Area_um2(i) = numel(cellPix) * pxArea_um2;

        % Meta
        Tn.Position(i) = s; Tc.Position(i) = s; Tw.Position(i) = s;
        Tn.NucleusID(i)= i; Tc.NucleusID(i)= i; Tw.NucleusID(i)= i;

        Tn.ThrChoice(i)= string(thrChoice); Tc.ThrChoice(i)= string(thrChoice); Tw.ThrChoice(i)= string(thrChoice);
        Tn.ThrInfo(i)  = string(thrInfo);   Tc.ThrInfo(i)  = string(thrInfo);   Tw.ThrInfo(i)  = string(thrInfo);
        Tn.CytSens(i)  = sensUsed;          Tc.CytSens(i)  = sensUsed;          Tw.CytSens(i)  = sensUsed;

        % ---- Save per-channel background value into CSVs (does NOT affect intensities) ----
        for c = 1:nC
            colName = sprintf('Bg_%s', matlab.lang.makeValidName(chanNames{c}));
            val = bgValue(c);
            if isnan(val)
                valNum = NaN;
            else
                valNum = val;
            end
            Tn.(colName)(i,1) = valNum;
            Tc.(colName)(i,1) = valNum;
            Tw.(colName)(i,1) = valNum;
        end

        % Intensities per channel
        for c = 1:nC
            nm = matlab.lang.makeValidName(chanNames{c});
            vN = channels{c}(nucPix);
            Tn.(sprintf('Mean_%s_nuc',nm))(i,1)   = mean(double(vN));
            Tn.(sprintf('IntDen_%s_nuc',nm))(i,1) = sum(double(vN));

            if ~isempty(cytoPix)
                vC = channels{c}(cytoPix);
                Tc.(sprintf('Mean_%s_cyto',nm))(i,1)   = mean(double(vC));
                Tc.(sprintf('IntDen_%s_cyto',nm))(i,1) = sum(double(vC));
            else
                Tc.(sprintf('Mean_%s_cyto',nm))(i,1)   = NaN;
                Tc.(sprintf('IntDen_%s_cyto',nm))(i,1) = NaN;
            end

            if ~isempty(cellPix)
                vW = channels{c}(cellPix);
                Tw.(sprintf('Mean_%s_cell',nm))(i,1)   = mean(double(vW));
                Tw.(sprintf('IntDen_%s_cell',nm))(i,1) = sum(double(vW));
            else
                Tw.(sprintf('Mean_%s_cell',nm))(i,1)   = NaN;
                Tw.(sprintf('IntDen_%s_cell',nm))(i,1) = NaN;
            end
        end
    end

    % Append tables
    T_nuc  = [T_nuc;  Tn]; %#ok<AGROW>
    T_cyto = [T_cyto; Tc]; %#ok<AGROW>
    T_cell = [T_cell; Tw]; %#ok<AGROW>

    % Masks
    nucRGB  = label2rgb(double(nucLabels),'jet','k','shuffle');
    cytRGB  = label2rgb(double(cytLabels),'jet','k','shuffle');
    cellRGB = label2rgb(double(cellLabels),'jet','k','shuffle');
    imwrite(nucRGB,  fullfile(outDir.Masks, sprintf('%s_pos%02d_nucleus.png',    opts.outputPrefix, s)));
    imwrite(cytRGB,  fullfile(outDir.Masks, sprintf('%s_pos%02d_cytoplasm.png',  opts.outputPrefix, s)));
    imwrite(cellRGB, fullfile(outDir.Masks, sprintf('%s_pos%02d_wholecell.png',  opts.outputPrefix, s)));

    % Segmentation overlay
    overlay = createSegOverlay(chNuc, nucLabels, cytLabels);
    imwrite(overlay, fullfile(outDir.Masks, sprintf('%s_pos%02d_overlay.png', opts.outputPrefix, s)));

    % Logging: per-position summary + per-channel background method/value
    logLines{end+1} = sprintf('Pos %02d | nuclei=%s (%s) | cyt_sens=%.3f | nucCount=%d', ...
        s, thrChoice, thrInfo, sensUsed, nCells); %#ok<AGROW>

    for c = 1:nC
        if isnan(bgValue(c))
            valStr = 'N/A';
        else
            valStr = sprintf('%.6f', bgValue(c));
        end
        logLines{end+1} = sprintf('    Channel %s: method=%s | bgValue=%s | reason=%s', ...
            chanNames{c}, bgMethod{c}, valStr, bgReason{c}); %#ok<AGROW>
    end
end

%% -------- Save CSVs (and create unified master) ----------
writetable(T_cell, fullfile(outDir.WholeCell, [opts.outputPrefix '_WholeCell.csv']));
writetable(T_nuc,  fullfile(outDir.Nucleus,   [opts.outputPrefix '_Nucleus.csv']));
writetable(T_cyto, fullfile(outDir.Cytoplasm, [opts.outputPrefix '_Cytoplasm.csv']));

U_cell = harmonizeForMaster(T_cell, 'WholeCell');
U_nuc  = harmonizeForMaster(T_nuc,  'Nucleus');
U_cyto = harmonizeForMaster(T_cyto, 'Cytoplasm');

masterALL = [U_cell; U_nuc; U_cyto];
writetable(masterALL, fullfile(outRoot, [opts.outputPrefix '_ALL.csv']));

for c = 1:nC
    nm = matlab.lang.makeValidName(chanNames{c});
    cols = {'Position','NucleusID','Area_um2','Compartment'};
    keepCols = {sprintf('Mean_%s',nm), sprintf('IntDen_%s',nm), 'ThrChoice','ThrInfo','CytSens'};

    tblN = U_nuc(:,  union(cols, keepCols));
    tblC = U_cyto(:, union(cols, keepCols));
    tblW = U_cell(:, union(cols, keepCols));

    tblN.Properties.VariableNames{'Area_um2'} = 'AreaNucleus_um2';
    tblC.Properties.VariableNames{'Area_um2'} = 'AreaCyto_um2';
    tblW.Properties.VariableNames{'Area_um2'} = 'AreaCell_um2';

    csvN = fullfile(outDir.Nucleus,   sprintf('%s_%s_Nucleus.csv',   opts.outputPrefix, nm));
    csvC = fullfile(outDir.Cytoplasm, sprintf('%s_%s_Cytoplasm.csv', opts.outputPrefix, nm));
    csvW = fullfile(outDir.WholeCell, sprintf('%s_%s_WholeCell.csv', opts.outputPrefix, nm));
    writetable(tblN, csvN); writetable(tblC, csvC); writetable(tblW, csvW);
end

%% -------- Write method log ----------
logPath = fullfile(outDir.Logs, 'method_log.txt');
fid = fopen(logPath, 'w');
if fid>0
    fprintf(fid, 'analyze_immunostain ? method log\n');
    fprintf(fid, 'Date: %s\n', datestr(now));
    fprintf(fid, 'PixelSize_um = %.6f (Area per pixel = %.8f um^2)\n', opts.pixelSize_um, pxArea_um2);
    fprintf(fid, 'CropRect = [%s]\n', num2str(cropRect));
    if isfield(bgMean,'value') && ~isempty(bgMean.value), fprintf(fid, 'Background (global scalar) = %.6f\n', bgMean.value); end
    if ~isempty(bgMean.ch405), fprintf(fid, 'Background ch405 = %.6f\n', bgMean.ch405); end
    if ~isempty(bgMean.ch488), fprintf(fid, 'Background ch488 = %.6f\n', bgMean.ch488); end
    if ~isempty(bgMean.ch647), fprintf(fid, 'Background ch647 = %.6f\n', bgMean.ch647); end
    fprintf(fid, 'Channels = %s\n', strjoin(chanNames, ', '));
    fprintf(fid, 'Options: useWatershed=%d, holeFill=%d, minNucleusArea_um2=%.3f\n', ...
        opts.useWatershed, opts.holeFill, opts.minNucleusArea_um2);
    fprintf(fid, '\nPer-position choices:\n');
    for i = 1:numel(logLines), fprintf(fid, '%s\n', logLines{i}); end
    fclose(fid);
end

%% -------- Cleanup ----------
try, reader.close(); catch, end
T = masterALL;
fprintf('Done. Results in: %s\n', outRoot);
end
% ======= end main =======

%% ================= helpers =================
function d = getDefaultOpts()
d.channelNames = struct('nuc405',1,'piezo488',2);
d.pixelSize_um = [];
d.gaussSigma = 1;
d.minNucleusArea_um2 = 20;
d.maxNucleusArea_um2 = inf;

d.askRenameChannels = true;
d.askOutputPrefix   = true;
d.askCrop           = true;
d.holeFill          = true;
d.useWatershed      = true;

d.cytoplasmSensitivity = 0.50;

d.background = struct('useROI',false,'value',[], 'ch405',[],'ch488',[],'ch647',[]);

d.outputPrefix = '';
d.saveBgSubChannels = true;  % keep default off; you can set true in opts
end

function L = splitAfterThreshold(nucMask, opts, pxArea_um2)
% Watershed ONLY for splitting touching nuclei (not full segmentation)
minPx = max(1, round(opts.minNucleusArea_um2/pxArea_um2));
nucMask = bwareaopen(nucMask, minPx);

% Distance transform
D = -bwdist(~nucMask);
D(~nucMask) = -Inf;

% Suppress shallow minima (avoid over-splitting)
D = imhmin(D, 1);

% Watershed
Lw = watershed(D);
L = Lw;
L(~nucMask) = 0;
end

function out = fillDefault(opts, defs)
out = defs;
f = fieldnames(opts);
for i = 1:numel(f), out.(f{i}) = opts.(f{i}); end
end

function [stack, px_um] = readMultichannel(reader, zPick)
if nargin<2 || isempty(zPick), zPick = 1; end
sizeX = reader.getSizeX();
sizeY = reader.getSizeY();
sizeZ = reader.getSizeZ();
sizeC = reader.getSizeC();
stack = zeros(sizeY, sizeX, sizeC, 'uint16');
tPick = 1;
zPick = min(max(1,zPick), sizeZ);
for c = 1:sizeC
    index = reader.getIndex(zPick-1, c-1, tPick-1);
    plane = bfGetPlane(reader, index+1);
    if ~isa(plane,'uint16'), plane = im2uint16(mat2gray(plane)); end
    stack(:,:,c) = plane;
end
px_um = [];
end

function satWarn(A, chanName)
if ~isfloat(A) && isa(A,'uint16')
    frac = nnz(A==65535)/numel(A);
    if frac > 0.01
        warning('~%.2f%%%% pixels saturated in channel %s', frac*100, string(chanName));
    end
end
end

function A = subtractBackground(A, bg)
if isempty(bg), return; end
A = double(A) - double(bg);
A(A<0) = 0;
A = cast(A, 'like', uint16(0));
A = uint16(A);
end

function bgMean = computeBackgroundMeans_fromImage(bgPath)
bgMean = struct('value',[],'ch405',[],'ch488',[],'ch647',[]);
try
    rdr = bfGetReader(bgPath);
    rdr.setSeries(0);
    [bgStack, ~] = readMultichannel(rdr, 1);
    try, rdr.close(); catch, end
    nbgC = size(bgStack,3);
    if nbgC >= 1, bgMean.ch405 = mean(double(bgStack(:,:,1)),'all'); end
    if nbgC >= 2, bgMean.ch488 = mean(double(bgStack(:,:,2)),'all'); end
    if nbgC >= 3, bgMean.ch647 = mean(double(bgStack(:,:,3)),'all'); end
catch
    try
        Ibg = imread(bgPath);
        m = mean(double(Ibg(:)));
        bgMean.ch405 = m; bgMean.ch488 = m; bgMean.ch647 = m;
    catch
        warning('Could not read background image: %s', bgPath);
    end
end
end

function J = normalize01(I)
I = double(I);
mn = min(I(:)); mx = max(I(:));
if mx>mn, J = (I-mn)/(mx-mn); else, J = zeros(size(I)); end
end

function [nucMask, nucLabels, choiceStr, infoStr] = segmentNucleiInteractive(I, opts, posIdx, pxArea_um2)
In = normalize01(I);
levOtsu = graythresh(In);

m1 = I > levOtsu * max(I(:));        % Otsu
m2 = I > (0.8*levOtsu) * max(I(:));  % Otsu relaxed
m3 = imbinarize(In, 'adaptive', 'Sensitivity', 0.5); % Adaptive

masks = {m1, m2, m3};
names = {'Otsu','Otsu x0.8','Adaptive(0.5)'};

minPx = max(1, round(opts.minNucleusArea_um2 / pxArea_um2));
for i = 1:numel(masks)
    if opts.holeFill, masks{i} = imfill(masks{i}, 'holes'); end
    masks{i} = bwareaopen(masks{i}, minPx);
end

% ---- LOCKED THRESHOLD METHOD (no popup) ----
nucMask   = masks{1};   % always Otsu
choiceStr = 'Otsu';
infoStr   = sprintf('Otsu=%.4f', levOtsu);

if isfinite(opts.maxNucleusArea_um2)
    CC = bwconncomp(nucMask);
    S  = regionprops(CC,'Area');
    keep = true(1,CC.NumObjects);
    for i=1:CC.NumObjects
        if (S(i).Area*pxArea_um2) > opts.maxNucleusArea_um2, keep(i)=false; end
    end
    L = labelmatrix(CC);
    nucMask = ismember(L, find(keep));
end

[nucLabels, ~] = bwlabel(nucMask);
end

function Lsplit = splitTouchingNuclei(nucMask, pxArea_um2, opts) %#ok<DEFNU>
minPx = max(1, round(opts.minNucleusArea_um2/pxArea_um2));
nucMask = bwareaopen(nucMask, minPx);
D = -bwdist(~nucMask);
D(~nucMask) = -Inf;
D = imhmin(D, 1);
L = watershed(D);
Lsplit = L; Lsplit(~nucMask) = 0;
end

function L = markerControlledWatershed(nucMask, opts, pxArea_um2)
minPx = max(1, round(opts.minNucleusArea_um2/pxArea_um2));
nucMask = bwareaopen(nucMask, minPx);
D = -bwdist(~nucMask);
D(~nucMask) = -Inf;
sureFG = imerode(nucMask, strel('disk',2));
markers = bwlabel(sureFG);
L = watershed(imimposemin(D, markers));
L(~nucMask) = 0;
end

function [cytMask, cytLabels, sensOut] = cytoUI_distanceConstrained(J, nucLabels, opts, posIdx) %#ok<INUSD>
% ---- LOCKED CYTOPLASM SENSITIVITY (always 0.8, no popup) ----
sensOut = 0.8;
[cytMask, cytLabels] = distanceExpansionCytoplasm(J, nucLabels, sensOut);
end

function [cytMask, cytLabels] = distanceExpansionCytoplasm(refChannel, nucLabels, sensitivity)
baseMask = imbinarize(normalize01(refChannel), 'adaptive', 'Sensitivity', sensitivity);
baseMask = imopen(baseMask, strel('disk',2));
baseMask = imclose(baseMask, strel('disk',2));
[~, idx] = bwdist(nucLabels > 0);
cytLabels = zeros(size(nucLabels), 'uint16');
cytLabels(baseMask) = uint16(nucLabels(idx(baseMask)));
cytLabels(nucLabels>0) = 0;
cytMask = cytLabels > 0;
end

function [stack, bgMethod, bgValue, bgReason] = enhancedBackgroundSubtraction(stack,bgMeans,chanNames)
% Strict background subtraction using ONLY the provided background image.
% No rolling-ball fallback is used here.

nC = size(stack,3);
bgMethod = cell(nC,1);
bgReason = cell(nC,1);
bgValue  = nan(nC,1);

% Compute a global fallback background if needed
allVals = [];
if ~isempty(bgMeans.ch405), allVals(end+1) = bgMeans.ch405; end
if ~isempty(bgMeans.ch488), allVals(end+1) = bgMeans.ch488; end
if ~isempty(bgMeans.ch647), allVals(end+1) = bgMeans.ch647; end
if isempty(bgMeans.value) && ~isempty(allVals)
    bgMeans.value = mean(allVals);  % Compute a global scalar from given channels
end

for c = 1:nC
    A = double(stack(:,:,c));
    nm = lower(chanNames{c});
    bg = [];
    method = '';
    reason = '';

    % 1) Wavelength-specific background match
    if contains(nm,'405') && ~isempty(bgMeans.ch405)
        bg = bgMeans.ch405;
        method = '405-specific background';
        reason = 'Channel name contains "405"';
    elseif contains(nm,'488') && ~isempty(bgMeans.ch488)
        bg = bgMeans.ch488;
        method = '488-specific background';
        reason = 'Channel name contains "488"';
    elseif contains(nm,'647') && ~isempty(bgMeans.ch647)
        bg = bgMeans.ch647;
        method = '647-specific background';
        reason = 'Channel name contains "647"';
    end

    % 2) Fallback to global background
    if isempty(bg) && ~isempty(bgMeans.value)
        bg = bgMeans.value;
        method = 'global scalar background';
        reason = 'No wavelength match';
    end

    % 3) If BG still empty ? no subtraction (bg=0)
    if isempty(bg)
        bg = 0;
        method = 'no subtraction';
        reason = 'Background image lacked usable channels';
    end

    % Apply subtraction
    A = A - double(bg);
    A(A < 0) = 0;
    stack(:,:,c) = uint16(A);

    % Save logs
    bgMethod{c} = method;
    bgReason{c} = reason;
    bgValue(c)  = bg;
end
end

function imwrite_uint16(path, I)
% Preserve full uint16 dynamic range exactly as in memory.
if ~isa(I,'uint16')
    I = uint16(I);   % No scaling, no mat2gray
end
imwrite(I, path, 'Compression','none');  % PNG/TIF support uint16
end

function sub = cropByRect(I, rect)
[H,W] = size(I);
x1 = max(1, rect(1)); y1 = max(1, rect(2));
w = max(1, rect(3));  h = max(1, rect(4));
x2 = min(W, x1 + w - 1);
y2 = min(H, y1 + h - 1);
sub = I(y1:y2, x1:x2);
end

function stackOut = cropStackByRect(stack, rect)
x1 = max(1, rect(1)); y1 = max(1, rect(2));
w = max(1, rect(3));  h = max(1, rect(4));
x2 = min(size(stack,2), x1 + w - 1);
y2 = min(size(stack,1), y1 + h - 1);
stackOut = stack(y1:y2, x1:x2, :);
end

function U = harmonizeForMaster(Tin, compartment)
U = Tin;
vars = U.Properties.VariableNames;
for k = 1:numel(vars)
    v = vars{k}; vnew = v;
    vnew = regexprep(vnew, '_nuc$',  '');
    vnew = regexprep(vnew, '_cyto$', '');
    vnew = regexprep(vnew, '_cell$', '');
    if ~strcmp(vnew, v)
        U.Properties.VariableNames{v} = vnew;
    end
end
if ~ismember('Compartment', U.Properties.VariableNames)
    U.Compartment = repmat(string(compartment), height(U), 1);
else
    U.Compartment(:) = string(compartment);
end
end

function RGB = createSegOverlay(base, nucLabels, cytLabels)
% Create an RGB overlay: grayscale base + yellow nucleus borders + cyan cyto borders.
base = mat2gray(base);
R = base; G = base; B = base;

perimN = bwperim(nucLabels > 0);   % nuclei borders
perimC = bwperim(cytLabels > 0);   % cyto borders

% Nuclei: yellow (R=1,G=1,B=0)
R(perimN) = 1; G(perimN) = 1; B(perimN) = 0;

% Cytoplasm: cyan (R=0,G=1,B=1)
R(perimC) = 0; G(perimC) = 1; B(perimC) = 1;

RGB = cat(3,R,G,B);
end
% 