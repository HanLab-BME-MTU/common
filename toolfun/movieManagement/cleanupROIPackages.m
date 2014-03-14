function cleanupROIPackages(MD, packageName, varargin)
%CLEANUPROIPACKAGES cleans packages in ROIs and optionally recreate new packages
%
%    cleanupROIPackages(MD, packageName) processes each ROIs from the
%    input movie tree and cleans a top-level package of type packageName.
%    It first finds a package of type packageName in the movie ancestor,
%    geos into each ROI and unlinks the package and all its processes from
%    the ROI if the package is shared.
%
%    cleanupROIPackages(MD, packageName, process2keep) additionally
%    preserves the processes of the top-level package of index specified by
%    process2keep from being unlinked. Additionally, it creates a new
%    package of type packageName in each ROI and links the kept process(es)
%    to each ROI package. Warning : this result in individual processes
%    being shared by multiple packages of different owners.
%
%    Examples:
%
%        cleanupROIPackages(MD, 'WindowingPackage')
%        cleanupROIPackages(MD, 'WindowingPackage', 1)
%
% Sebastien Besson, Mar 2014

% Input check
ip = inputParser();
ip.addRequired('packageName', @ischar);
ip.addOptional('process2keep', [], @isnumeric);
ip.parse(packageName, varargin{:})

% Retrieve package of input class from ancestor
ancestor = MD.getAncestor();
packageIndex = ancestor.getPackageIndex(packageName, 1, false);
package = ancestor.getPackage(packageIndex);

% Unlink package and processes from children movies
for movie  = ancestor.getDescendants()
    cleanPackage(movie, package, ip.Results.process2keep)
end

% Return if no process is kept
if isempty(ip.Results.process2keep), return; end

% Recreate packages and link kept process to each of them
for movie  = ancestor.getDescendants()
    recreatePackage(movie, package, ip.Results.process2keep)
end

function cleanPackage(movie, package, process2keep)

% Compute list of processes to unlink
processes2clean = ~cellfun(@isempty, package.processes_);
processes2clean(process2keep) = false;

% Unlink processes and package
for i = find(processes2clean)
    movie.unlinkProcess(package.getProcess(i));
end
movie.unlinkPackage(package);

function recreatePackage(movie, package, process2keep)

% Create a new package using the default constructor
packageConstr = str2func(class(package));
newPackage = packageConstr(movie);
movie.addPackage(newPackage);

% Link parent processes
for i = find(process2keep)
    newPackage.setProcess(i, package.getProcess(i));
end