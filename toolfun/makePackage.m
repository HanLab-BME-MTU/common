function makePackage(outDir)
%MAKEPACKAGE builds the list of selected packages
% 
% makePackage(outDir)
% 
% This function copies all the files needed to run the selected packages data
% processing package into a single folder for upload to the website. It
% assumes you have checked out all the files from the SVN repository, and
% that they are all in your matlab path.
% 
% Input:
% 
%   outDir - The directory to copy all the package files to.
%
%
% Sebastien Besson
% 03/2011

% List of all available packages
fullPackageList={'SegmentationPackage';'BiosensorsPackage';'UTrackPackage'};

% Ask the user which packages to include
[packageIndx,ok] = listdlg('PromptString','Select the package(s) to compile:',...
    'SelectionMode','multiple','ListString',fullPackageList);
if ~ok, return; end
packageList = fullPackageList(packageIndx);

% Retrieve the full path of the corresponding packages
packageLocation=cellfun(@which,packageList,'UniformOutput',false);
if any(cellfun(@isempty,packageLocation)),
    errordlg('Could not locate on of the selected packages.')
    return;
end
packageDir=cellfun(@fileparts,packageLocation,'UniformOutput',false);
 
% Ask for the output directory if not supplied
if nargin < 1 || isempty(outDir)    
    outDir = uigetdir(pwd,'Select output dir:');
end
    
%Get m files from packages directory
disp('Getting file list...')
packageFuns= cellfun(@(x) dir([x filesep '*.m']),packageDir,'UniformOutput',false);
packageFuns=vertcat(packageFuns{:});
packageFuns = {packageFuns(:).name}';
%Remove this function from the list if present (but it shouldn't!)
packageFuns(strcmp(packageFuns,mfilename)) = [];

%Add the few odd files that are in other areas of common
if any(strcmp(packageList,'UTrackPackage'))
    packageFuns = vertcat(packageFuns,{'kalmanResMemLM.m';'kalmanInitLinearMotion.m';...
        'kalmanGainLinearMotion';'kalmanReverseLinearMotion';'costMatLinearMotionLink2';'costMatLinearMotionCloseGaps2'});
end

if any(strcmp(packageList,'BiosensorsPackage'))
    packageFuns = vertcat(packageFuns,{'refineMovieMasks.m','separateNumberedFiles.m'}');
end

%Get everything these functions depend on also
packageFuns = depfun_notoolbox(packageFuns);

% Find associated documentation files
[packageFunsPaths packageFunsNames]=cellfun(@fileparts,packageFuns,'UniformOutput',false);
isDocFile = logical(cellfun(@(x,y) exist([x filesep 'doc' filesep y '.pdf'],'file'),...
    packageFunsPaths,packageFunsNames));
packageDocFiles = cellfun(@(x,y) [x filesep 'doc' filesep y '.pdf'],...
    packageFunsPaths(isDocFile),packageFunsNames(isDocFile),'UniformOutput',false);

%Check and display toolbox dependency
disp('Checking toolbox dependency...')
tbs = toolboxesUsed(packageFuns);
disp('The package uses the following toolboxes:')
disp(tbs)

% Check for the presence of GUI fig files and append them
isGUIFile =logical(cellfun(@(x) exist([x(1:end-2) '.fig'],'file'),packageFuns));
packageFigs = cellfun(@(x) [x(1:end-2) '.fig'],packageFuns(isGUIFile),'UniformOutput',false);

packageIcons =  which('lccbGuiIcons.mat');
% Add icon files
packageFiles=vertcat(packageFuns,packageFigs,packageIcons);

% Create package output directory if non-existing
if ~exist(outDir,'dir'), mkdir(outDir); end

nFiles = numel(packageFiles);
disp(['Copying all '  num2str(nFiles) ' files ...'])
for j = 1:nFiles
    iLFS = max(regexp(packageFiles{j},filesep));
    copyfile(packageFiles{j},[outDir filesep packageFiles{j}(iLFS+1:end)]);
end

% Create doc output directory if non-existing
docDir=[outDir filesep 'doc'];
if ~exist(docDir,'dir'), mkdir(docDir); end

nDocFiles = numel(packageDocFiles);
disp(['Copying all '  num2str(nDocFiles) ' files ...'])

for nfile=1:numel(nDocFiles)
    iLFS = max(regexp(packageDocFiles{nfile},filesep));
    copyfile(packageDocFiles{nfile},[docDir filesep packageDocFiles{nfile}(iLFS+1:end)]);
end
    
disp(['Finished. Wrote package to ' outDir])
