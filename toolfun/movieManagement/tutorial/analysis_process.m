%% Analysis set-up via Command line


%% Initialization

dvPath = '~/Desktop/2014Mar26/Actin-TM2.ome.tif';
MD = MovieData.load(dvPath);
fprintf(1, 'Object saved under: %s\n', MD.getFullPath());
fprintf(1, 'Output directory for analysis: %s\n', MD.outputDirectory_);

%% Reset

% Reset analysis
MD.reset();

%% First process

% Set-up process via command line interace
process = ThresholdProcess(MD);
MD.addProcess(process);
processIndex = process.getIndex();
parameters = process.getParameters();

fprintf(1, 'Process %g: %s\n', processIndex, process.getName());

% Retrieve analysis parameters
disp('Parameters');
disp(parameters);
disp('  Default parameters:');
disp(MD.getProcess(processIndex).getDefaultParams(MD));

% Initial analysis status
disp('Status');
fprintf(1,'  Process has been run successfully: ');
if process.success_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Parameters have been modified since last successful run: ');
if process.procChanged_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Input has been updated by an upstream process: ');
if ~process.updated_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end

% Run the first process
process.run();

% Post-run status
disp('Status');
fprintf(1,'  Process has been run successfully: ');
if process.success_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Parameters have been modified since last successful run: ');
if process.procChanged_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Input has been updated by an upstream process: ');
if ~process.updated_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end

%% Parameters modification

parameters = process.getParameters();

%
parameters.GaussFilterSigma = 1;
disp('Setting new parameters');
process.setParameters(parameters);

% Post parameters modification analysis status
disp('Status');
fprintf(1,'  Process has been run successfully: ');
if process.success_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Parameters have been modified since last successful run: ');
if process.procChanged_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Input has been updated by an upstream process: ');
if ~process.updated_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end

% Second run of the first process
process.run();

% Post parameters modification analysis status
disp('Status');
fprintf(1,'  Process has been run successfully: ');
if process.success_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Parameters have been modified since last successful run: ');
if process.procChanged_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end
fprintf(1,'  Input has been updated by an upstream process: ');
if ~process.updated_, fprintf(1, 'yes\n'); else fprintf(1, 'no\n'); end

