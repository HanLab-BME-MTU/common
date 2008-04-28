function [] = runLifetimeAnalysis(restrict,shape)

% runLifetimeAnalysis takes all data folders under a given condition and fills in
% missing tracking and lifetime data
%
% SYNOPSIS [] = runLifetimeAnalysis()
%
% INPUT     restrict    = (optional) time restriction in seconds, e.g. 300
%           shape   = (optional) shape for populations, e.g. [2 2 1], where 1
%                   indicates an exponential distribution, and 2 a Rayleigh
%                   distribution
%
% OUTPUT
%
% REMARKS   When specifying either the restrict or the shape parameters and
%           not specifying the other parameter leave the non specified
%           parameter as an empty field (e.g. runLiftimeAnalysis(300,[]) or
%           runLifetimeAnalysis([],[2 2 1])
%
% Daniel Nunez, March 5, 2008

oldDir = cd;

%Load Data
%ask user to specify movies to analyse
[experiment] = loadIndividualMovies();
%go to condition folder
cd(experiment(1).source);
cd ..
cd ..
%ask user to identify folder to save data to; condition folder is used as
%default
directory = uigetdir(pwd,['Specify folder to store lifetime analysis result.\n If empty default (' cd ') will be used']);
%get directory name to use as default for analysis results
dirName = findstr(directory,filesep);
dirName = directory(dirName(end)+1:end);
%if user doesn't input anything, use default
if ~ischar(directory)
    directory = pwd;
end
%Ask user to name file
fileName = input('Specify name for lifetime analysis result files (date will be included automatically).','s');
if isempty(fileName)
    fileName = dirName;
end


%Default Values
%if user did not specify inputs then use defaults
if nargin == 0
    restrict = 300; %restrict pit lifetimes
    shape = [2 2 1]; %distribution shape vector
    %if restrict is specified, but empty then use default
elseif isempty(restrict) && length(shape) == 3
    restrict = 300;
    %if shape is specified, but empty then use default
elseif isempty(shape) && length(restrict) == 1
    shape = [2 2 1];
    %if two inputs do nothing unless restrict is a vector and shape a scalar
elseif nargin == 2
    if length(shape) == 1 && length(restrict) == 3
        new_shape = restrict;
        new_restrict = shape;
        shape = new_shape;
        restrict = new_restrict;
        warning('Inputs may have been misinterpreted.');
    end
else
    error('Input must of the form (),(restrict,[]), ([],shape), or (restrict,shape).\n Restrict must be of length 1 and shape of length 3.')
end

%Fill in movie length and image size
[experiment] = determineMovieLength(experiment);
[experiment] = determineImagesize(experiment);
%Create lftInfo file if not already there
[experiment] = determineLifetimeInfo(experiment);
%load lftInfo into experiment structures
for iexp = 1:length(experiment)
    cd([experiment(iexp).source filesep 'LifetimeInfo']);
    load('lftInfo.mat');
    experiment(iexp).lftInfo = lftInfo;
end

%Analyse Data
[compactRes, data] = lifetimeCompactFitData(experiment, restrict, shape);

%Save Data
secureSave([directory filesep fileName datestr(now,'yyyymmdd')],'compactRes');

fid = fopen([directory filesep fileName datestr(now,'yyyymmdd') '.txt'],'wt');
fprintf(fid,'%s\n',[fileName ',shape,C,deltaC,Tau,deltaTau,Tau50']);
fprintf(fid,'%s\n',['PP,' num2str(0) ',' num2str(compactRes.contr(1)) ',' num2str(compactRes.contrError(1)) ',' num2str(compactRes.tau(1)) ',' num2str(compactRes.tauError(1)) ',' num2str(compactRes.tau50(1)) ',' num2str(compactRes.numcells) ',cell number']);
fprintf(fid,'%s\n',['P1,' num2str(shape(1)) ',' num2str(compactRes.contr(2)) ',' num2str(compactRes.contrError(2)) ',' num2str(compactRes.tau(2)) ',' num2str(compactRes.tauError(2)) ',' num2str(compactRes.tau50(2)) ',' num2str(compactRes.numtraj) ',traj number']);
fprintf(fid,'%s\n',['P2,' num2str(shape(2)) ',' num2str(compactRes.contr(3)) ',' num2str(compactRes.contrError(3)) ',' num2str(compactRes.tau(3)) ',' num2str(compactRes.tauError(3)) ',' num2str(compactRes.tau50(3)) ',' 'num2str(compactRes.numcells)' ',density']);
fprintf(fid,'%s\n',['P3,' num2str(shape(3)) ',' num2str(compactRes.contr(4)) ',' num2str(compactRes.contrError(4)) ',' num2str(compactRes.tau(4)) ',' num2str(compactRes.tauError(4)) ',' num2str(compactRes.tau50(4)) ',' 'num2str(compactRes.numcells)' ',delta density']);
fprintf(fid,'%s\n',['Max lifetime for analysis (restrict):,' num2str(restrict)]);
fprintf(fid,'%s\n',['Movies Used:']);
for iexperiment = 1:length(experiment)
fprintf(fid,'%s\n',[experiment(iexperiment).source]);
end
fclose(fid);

cd(oldDir);

end %of function