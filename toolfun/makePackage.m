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
%
% Sebastien Besson, July 2011

% List available packages
packages(5)=struct('name','','initFiles','','additionalFiles','');
packages(1).name='Segmentation';
packages(1).initFiles={'SegmentationPackage'};
packages(2).name='Biosensors';
packages(2).initFiles={'BiosensorsPackage'};
packages(3).name='U-Track';
packages(3).initFiles={'UTrackPackage'};
packages(4).name='QFSM';
packages(4).initFiles={'QFSMPackage'};
for i=1:4, packages(i).additionalFiles = {'lccbGUIicons.mat'}; end
packages(5).name='plusTipTracker';
packages(5).initFiles={'plusTipGetTracks';'plusTipSeeTracks';...
    'plusTipParamSweepGUI';'plusTipGroupAnalysis'};
packages(5).additionalFiles = {'pTT_logo_sm.png';'help_icon.png'};

% Retrieve the list of valid packages
isValidPackage=@(x)all(cellfun(@(y)exist(y,'file'),packages(x).initFiles));
validPackages=arrayfun(isValidPackage,1:numel(packages));
packages(~validPackages)=[];
if isempty(packages), error('No package found'); end

% Ask the user which packages to build
[packageIndx,status] = listdlg('PromptString','Select the package(s) to build:',...
    'SelectionMode','multiple','ListString',{packages(:).name});
if ~status, return; end
packages=packages(packageIndx);

% Ask for the output directory if not supplied
if nargin < 1 || isempty(outDir)    
    outDir = uigetdir(pwd,'Select output dir:');
end
    
%Get all the function dependencies and display toolboxes
[packageFuns toolboxesUsed] = getFunDependencies(vertcat(packages.initFiles));
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

% List all files which are neither M-files nor FIG-files nor MAT-files
uniquePackageFunsExt = unique(packageFunsExt);
matExt={'.fig';'.m';'.mat'};
matExtIndx = ismember(uniquePackageFunsExt,matExt);
mexExt=uniquePackageFunsExt(~matExtIndx);
mexFunsIndx = find(ismember(packageFunsExt,mexExt));

% Get all files in the same folder as these found MEX-files
packageMexList=arrayfun(@(x)  dir([packageFunsPaths{x} filesep '*.*']),...
    mexFunsIndx,'Unif',false);
packageMexFunsPaths=packageFunsPaths(mexFunsIndx);
packageMexFunsNames = @(x) strcat([packageMexFunsPaths{x} filesep],...
    {packageMexList{x}(~[packageMexList{x}.isdir]).name}');
packageMexFuns = arrayfun(@(x) packageMexFunsNames(x),1:numel(mexFunsIndx),'Unif',false);
packageMexFuns =vertcat(packageMexFuns{:});

% Remove C-files
if ~isempty(packageMexFuns)
    cFiles=~cellfun(@isempty,regexp(packageMexFuns,'.c$','once'));
    packageMexFuns(cFiles)=[];
end

% Add additional files (not included in dependencies)
packageIcons = unique(cellfun(@which,vertcat(packages.additionalFiles),'UniformOutput',false));

% Concatenate all matlab files but the documentation
packageFiles=vertcat(packageFuns,packageFigs,packageIcons,packageMexFuns);

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
disp('Creating/cleaning release documentation directory...')
docDir=[outDir filesep 'doc'];
mkClrDir(docDir);

% Copy documentation files
nDocFiles = numel(packageDocs);
disp(['Copying all '  num2str(nDocFiles) ' files ...'])
for nfile=1:numel(packageDocs)
    iLFS = max(regexp(packageDocs{nfile},filesep));
    copyfile(packageDocs{nfile},[docDir filesep packageDocs{nfile}(iLFS+1:end)]);
end
    
disp(['Finished. Wrote package to ' outDir])
