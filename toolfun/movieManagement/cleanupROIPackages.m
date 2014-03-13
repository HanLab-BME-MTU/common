function cleanupROIPackages(MD, packageName, varargin)

% Input check
ip = inputParser();
ip.addRequired('packageName', @ischar);
ip.addOptional('process2keep', [], @isnumeric);
ip.parse(packageName, varargin{:})

% Retrieve package of input class from ancestor
ancestor = MD.getAncestor();
packageIndex = ancestor.getPackageIndex('WindowingPackage', 1, false);
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