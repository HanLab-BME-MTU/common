function makePackage(outDir)
% Build the list of selected packages
% 
% makePackage(outDir)
% 
% This function copies all the files needed to run the selected packages data
% processing package into a single folder for upload to the website. It
% assumes you have checked out all the files from the SVN repository, and
% that they are all in your matlab path.
% 
% INPUT:
% 
%   outDir - The directory to copy all the package files to.
%

% Sebastien Besson, July 2011 (last modified: Jul 2012)

% List available packages and additional files required for running them
packageList = {
    'BiosensorsPackage';...
    'IntegratorPackage'
    'QFSMPackage'
    'SegmentationPackage'
    'TrackingPackage'
    'WindowingPackage'};
validPackage = cellfun(@(x) exist(x,'class')==8,packageList);
packageList = packageList(validPackage);
if isempty(packageList), error('No package found'); end

% Create package names
packageNames = cellfun(@(x) eval([x '.getName']),packageList,'Unif',0);
[packageNames,index]=sort(packageNames);
packageList=packageList(index);

% Ask the user which packages to build
[packageIndx,status] = listdlg('PromptString','Select the package(s) to build:',...
    'SelectionMode','multiple','ListString',packageNames);
if ~status, return; end
packageList=packageList(packageIndx);

% Ask for the output directory if not supplied
if nargin < 1 || isempty(outDir)    
    outDir = uigetdir(pwd,'Select output dir:');
end
    
%Get all the function dependencies and display toolboxes
[packageFuns toolboxesUsed] = getFunDependencies(packageList);
disp('The package uses the following toolboxes:')
disp(toolboxesUsed)

%% Additional files can be found under four types of format:
%   * GUIs may have associated *.fig
%   * Processes and GUIs may have associated *.pdf
%   * Mex-files have many extensions depending on the OS
%   * Icons or other MAT-files

% Split functions into paths, filenames and extensions for search
[packageFunsPaths packageFunsNames packageFunsExt]=...
    cellfun(@fileparts,packageFuns,'UniformOutput',false);

% Find associated documentation files
isDocFile = logical(cellfun(@(x,y) exist([x filesep 'doc' filesep y '.pdf'],'file'),...
    packageFunsPaths,packageFunsNames));
packageDocs = cellfun(@(x,y) [x filesep 'doc' filesep y '.pdf'],...
    packageFunsPaths(isDocFile),packageFunsNames(isDocFile),'UniformOutput',false);

% Get GUI fig files
isGUIFile =logical(cellfun(@(x) exist([x(1:end-2) '.fig'],'file'),packageFuns));
packageFigs = cellfun(@(x) [x(1:end-2) '.fig'],packageFuns(isGUIFile),'UniformOutput',false);

% List all mex files in the package
mexFunsExt={'.dll';'.mexglx';'.mexmaci';'.mexmaci64';'.mexa64';'.mexw64'};
mexFunsIndx= find(ismember(packageFunsExt,mexFunsExt));
packageMexList=arrayfun(@(x)  dir([packageFunsPaths{x} filesep packageFunsNames{x} '.*']),...
    mexFunsIndx,'Unif',false);
packageMexFunsPaths=packageFunsPaths(mexFunsIndx);
packageMexFunsNames = @(x) strcat([packageMexFunsPaths{x} filesep],...
    {packageMexList{x}(~[packageMexList{x}.isdir]).name}');
packageMexFuns = arrayfun(@(x) packageMexFunsNames(x),1:numel(mexFunsIndx),'Unif',false);
packageMexFuns =vertcat(packageMexFuns{:});

% Remove additional compilation files
if ~isempty(packageMexFuns)
    compFunsExt = {'.c';'.cpp';'.h';'.nb';'.m'};
    for i=1:numel(compFunsExt)
        indx = ~cellfun(@isempty,regexp(packageMexFuns,[compFunsExt{i} '$'],'once'));
        packageMexFuns(indx)=[];
    end
end

% Get the main path to the icons folder
iconsPath = fullfile(fileparts(which('packageGUI.m')),'icons');
icons = dir([iconsPath filesep '*.png']);
packageIcons = arrayfun(@(x) [iconsPath filesep x.name],icons,'Unif',false);

% Concatenate all matlab files but the documentation
packageFiles=vertcat(packageFuns,packageFigs,packageMexFuns);

%% Export package files
% Create package output directory if non-existing
disp('Creating/cleaning release directory...')
mkClrDir(outDir);

% Copy function files
nFiles = numel(packageFiles);
disp(['Copying all '  num2str(nFiles) ' files ...'])
for j = 1:nFiles
    iLFS = max(regexp(packageFiles{j},filesep));
    copyfile(packageFiles{j},[outDir filesep packageFiles{j}(iLFS+1:end)]);
end

% Create doc output directory if non-existing
disp('Creating/cleaning release icons directory...')
iconsDir=[outDir filesep 'icons'];
mkClrDir(iconsDir);

% Copy documentation files
nIcons = numel(packageIcons);
disp(['Copying all '  num2str(nIcons) ' files ...'])
for i = 1:numel(packageIcons)
    iLFS = max(regexp(packageIcons{i},filesep));
    copyfile(packageIcons{i}, [iconsDir filesep packageIcons{i}(iLFS+1:end)]);
end

% Create doc output directory if non-existing
disp('Creating/cleaning release documentation directory...')
docDir=[outDir filesep 'doc'];
mkClrDir(docDir);

% Copy documentation files
nDocFiles = numel(packageDocs);
disp(['Copying all '  num2str(nDocFiles) ' files ...'])
for i = 1 : numel(packageDocs)
    iLFS = max(regexp(packageDocs{i},filesep));
    copyfile(packageDocs{i},[docDir filesep packageDocs{i}(iLFS+1:end)]);
end
    
disp(['Finished. Wrote package to ' outDir])
