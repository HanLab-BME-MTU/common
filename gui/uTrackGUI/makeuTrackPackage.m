function makeuTrackPackage(outDir,packageDir,comDir)
%MAKEUTRACKPACKAGE builds the uTrack package
% 
% makeuTrackPackage(outDir,uTrackDir,comDir)
% 
% THis function copies all the files needed to run the uTrack data
% processing package into a single folder for upload to the website. It
% assumes you have checked out all the files from the SVN repository, and
% that they are all in your matlab path.
% 
% Input:
% 
%   outDir - The directory to copy all the package files to.
%
%   packageDir - the directory you have checked out the package
% 
%   comDir - The directory you have checked out the common project to.
% 
% Hunter Elliott
% 12/2010
%
% Modified by Sebastien Besson
% 03/2011

if nargin < 1 || isempty(outDir)    
    outDir = uigetdir(pwd,'Select output dir:');
end
if nargin < 2 || isempty(packageDir)    
    packageDir = uigetdir(pwd,'Select package dir:');
end
if nargin < 3 || isempty(comDir)    
    comDir = uigetdir(pwd,'Select common dir:');
end

    
disp('Getting file list...')

%Get m files from package directory
packageFuns={'/home/sb286/Documents/MATLAB/common/gui/uTrackGUI/uTrackPackageGUI.m'};
% packageFuns = dir([packageDir filesep '*.m']);
% packageFuns = {packageFuns(:).name}';
% %Remove this function from the list. Fucking recursion.
% packageFuns(strcmp(packageFuns,mfilename)) = [];


%Get everything these depend on also
packageFuns = depfun_notoolbox(packageFuns);
%Get extra non-function files (figures for GUIs)
packageExtras = dir([packageDir filesep '*.fig']);
packageExtras = {packageExtras(:).name}';
packageExtras = cellfun(@(x)([packageDir filesep x]),packageExtras,'UniformOutput',false);

%Get the movie management functions
comDir = [comDir filesep 'toolfun' filesep 'movieManagement'];
comFuns = dir([comDir filesep '*.m']);
comFuns = {comFuns(:).name}';
%Add the few odd files that are in other areas of common
comFuns = vertcat(comFuns,{'kalmanResMemLM.m';'kalmanInitLinearMotion.m';...
    'kalmanGainLinearMotion';'kalmanReverseLinearMotion';'costMatLinearMotionLink2';'costMatLinearMotionCloseGaps2'});

%..and everything these functions depend on
comFuns = depfun_notoolbox(comFuns);

%Get fig files from common also
comExtras = dir([comDir filesep '*.fig']);
comExtras = {comExtras(:).name}';
comExtras{end+1} = 'lccbGuiIcons.mat';
comExtras = cellfun(@(x)([comDir filesep x]),comExtras,'UniformOutput',false);

%Check and display toolbox dependency
disp('Checking toolbox dependency...')
tbs = toolboxesUsed(vertcat(packageFuns,comFuns));

disp('The package uses the following toolboxes:')
disp(tbs)

if ~exist(outDir,'dir')
    mkdir(outDir)
end

allFiles = vertcat(packageFuns,packageExtras,comFuns,comExtras);
nFiles = numel(allFiles);

disp(['Copying all '  num2str(nFiles) ' files ...'])

for j = 1:nFiles

    iLFS = max(regexp(allFiles{j},filesep));
    copyfile(allFiles{j},[outDir filesep allFiles{j}(iLFS+1:end)]);

end

% Get documentation files and copy them to a doc directory (still very crude)
packageFunsDirs=unique(cellfun(@getFilenameBody,packageFuns,'UniformOutput',false));
isFunsDocDir=cellfun(@(x) exist([x filesep 'doc'],'dir')==7,packageFunsDirs);

docDir=[outDir filesep 'doc'];
if ~exist(docDir,'dir')
    mkdir(docDir)
end

allDocFiles={};
packageDocsDirs=cellfun(@(x) [x filesep 'doc'],{packageFunsDirs{isFunsDocDir}},'UniformOutput',false);
for j = 1:length(packageDocsDirs)
	docList=dir([packageDocsDirs{j} filesep '*.pdf']);
    docFiles=arrayfun(@(x)([packageDocsDirs{j} filesep x.name]),docList,'UniformOutput',false);
    allDocFiles=vertcat(allDocFiles,docFiles);
end

for nfile=1:length(allDocFiles)
    iLFS = max(regexp(allDocFiles{nfile},filesep));
    copyfile(allDocFiles{nfile},[docDir filesep allDocFiles{nfile}(iLFS+1:end)]);
end
    
disp(['Finished. Wrote package to ' outDir])
