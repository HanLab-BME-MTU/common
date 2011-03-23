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
% $LastChangedBy::                  $
% $LastChangedDate::                 $

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
    packageFuns = vertcat(packageFuns,{});
end

if any(strcmp(packageList,'BiosensorsPackage'))
    packageFuns = vertcat(packageFuns,{'refineMovieMasks.m','separateNumberedFiles.m'}');
end

%Get everything these functions depend on also
packageFuns = depfun_notoolbox(packageFuns);

%Check and display toolbox dependency
disp('Checking toolbox dependency...')
tbs = toolboxesUsed(packageFuns);
disp('The package uses the following toolboxes:')
disp(tbs)

% Additional functions can be found under three types of format:
%   * GUIs may have associated *.fig
%   * Processes and GUIs may have associated *.pdf
%   * Mex-files can be found under many format depending on the OS
% Split functions into paths, filenames and extensions for search
[packageFunsPaths packageFunsNames packageFunsExt]=...
    cellfun(@fileparts,packageFuns,'UniformOutput',false);


% Find associated documentation files
isDocFile = logical(cellfun(@(x,y) exist([x filesep 'doc' filesep y '.pdf'],'file'),...
    packageFunsPaths,packageFunsNames));
packageDocs = cellfun(@(x,y) [x filesep 'doc' filesep y '.pdf'],...
    packageFunsPaths(isDocFile),packageFunsNames(isDocFile),'UniformOutput',false);

% Get GUI-associated fig files
isGUIFile =logical(cellfun(@(x) exist([x(1:end-2) '.fig'],'file'),packageFuns));
packageFigs = cellfun(@(x) [x(1:end-2) '.fig'],packageFuns(isGUIFile),'UniformOutput',false);


% Add all versions of MEX files
uniquePackageFunsExt = unique(packageFunsExt);
matExt={'.fig';'.m';'.mat'};
matExtIndx = ismember(uniquePackageFunsExt,matExt);
mexExt=uniquePackageFunsExt(~matExtIndx);
mexFunsIndx = find(ismember(packageFunsExt,mexExt));

% List all files in the same folder as found MEX-files
packageMexList=arrayfun(@(x)  dir([packageFunsPaths{x} filesep '*.*']),...
    mexFunsIndx,'Unif',false);
packageMexFunsPaths=packageFunsPaths(mexFunsIndx);
packageMexFunsNames = @(x) strcat([packageMexFunsPaths{x} filesep],...
    {packageMexList{x}(~[packageMexList{x}.isdir]).name}');
packageMexFuns = arrayfun(@(x) packageMexFunsNames(x),1:numel(mexFunsIndx),'Unif',false);
packageMexFuns =vertcat(packageMexFuns{:});


% Add icons
packageIcons =  which('lccbGuiIcons.mat');

% Concatenate all files but the documentation
packageFiles=vertcat(packageFuns,packageFigs,packageIcons,packageMexFuns);

% Create package output directory if non-existing
disp('Creating/cleaning release directory...')
mkClrDir(outDir);

nFiles = numel(packageFiles);
disp(['Copying all '  num2str(nFiles) ' files ...'])
for j = 1:nFiles
    iLFS = max(regexp(packageFiles{j},filesep));
    copyfile(packageFiles{j},[outDir filesep packageFiles{j}(iLFS+1:end)]);
end

% Create doc output directory if non-existing
disp('Creating/cleaning release directory...')
docDir=[outDir filesep 'doc'];
mkClrDir(docDir);

nDocFiles = numel(packageDocs);
disp(['Copying all '  num2str(nDocFiles) ' files ...'])

for nfile=1:numel(packageDocs)
    iLFS = max(regexp(packageDocs{nfile},filesep));
    copyfile(packageDocs{nfile},[docDir filesep packageDocs{nfile}(iLFS+1:end)]);
end
    
disp(['Finished. Wrote package to ' outDir])
