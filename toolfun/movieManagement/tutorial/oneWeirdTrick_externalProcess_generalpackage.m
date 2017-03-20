% This example is meant to demonstrate the use of ExternalProcess and
% GenericPackage. These two classes are meant to allow the developer to
% quickly take advantage of the Process and Package APIs while prototyping
% an arbitrary function. In particular, with just a few commands, the
% Package GUI routines are made accessible.
%
% Below we combine a fully developed Process, ThresholdProcess, with a
% prototype process that says "Hello world!" into a GenericPackage. We then
% show the Package GUI where the the Processes can be configured and run.
%
% Author(s)
% Andrew Jamieson
% Mark Kittisopikul
% March 2017

%% Configuration
config.imageData = 'example.tif';

%% Load arbitrary image dataset
MD = MovieData.load(config.imageData);

%% Create and configure processes
% ThresholdProcess is a predefined, developed process with default
% parameters
threshProc = ThresholdProcess(MD);

% ExternalProcess is used to prototype our Hello World function
extProc = ExternalProcess(MD,'Say something',@(p) disp(p.getParameters().text));
% ExternalProcess uses a struct with no fields as a default
% Thus, we need to set the parameters with the field text
extProc.setParameters(struct('text','Hello world!'));

%% Add Processes and setup to MovieData object
MD.addProcess(threshProc);
MD.addProcess(extProc);
% This creates a GenericPackage with all the Processes contained in MD
% However, you can also give an arbitrary set of Processes
MD.addPackage(GenericPackage(MD));

% You may also want to configure the dependency matrix
% openvar('MD.packages_{1}.dependencyMatrix_')
% Or asssign a name
% MD.packages_{1}.name_ = 'Hello World';


%% Show GUI
MD.packages_{1}.GUI(MD);
% Alternatively, and with support for multiple GenericPackges use
% MD.packages_{1}.showGUI()
