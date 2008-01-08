function duplicateList = checkDuplicateFiles(matlabHome)
%CHECKDUPLICATEFILES checks for duplicate files in MATLAB-HOME
%
% SYNOPSIS: duplicateList = checkDuplicateFiles(matlabHome)
%
% INPUT matlabHome (optional): Directory where all the matlab files are stored
%
% OUTPUT duplicateList: n-by-2 cell array with {fileName, pathName}
%
% REMARKS This is a bit of a hack
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 30-Nov-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECK INPUT

duplicateList = cell(0,2);

if nargin == 0 || isempty(matlabHome)
    % check for MATLABHOME, HOME
    matlabHome = getenv('MATLABHOME');
    if isempty(matlabHome)
        % lccb-compatibility
        matlabHome = getenv('HOME');
    end
    if isempty(matlabHome)
        % warn
        warning('(MATLAB)HOME is not defined. Supply path')
        return
    end
    matlabHome = fullfile(matlabHome,'matlab');
end

%% find duplicate m-files

% get files
lof = searchFiles('\.m$','',matlabHome);

% remove @, remove contents.m
for i = length(lof):-1:1,
    if any(findstr(lof{i,2},'@')) || strcmpi(lof{i,1},'contents.m'),
        lof(i,:)=[];
    end,
end

% find unique filenames
names=strvcat(lof{:,1});
[un,m,n]=unique(names,'rows');
[ue,num]=countEntries(n);

% check for duplicates
idx=find(num>1);
for ii=idx',
    duplicateList = [duplicateList;lof(n==ue(ii),:)];
end