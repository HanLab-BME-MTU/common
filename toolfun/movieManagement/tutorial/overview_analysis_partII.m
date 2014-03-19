%% Analysis set-up via Command line
omeTiffPath = '~/Desktop/2014Mar20/case1_higherSNR.ome.tiff';
MD = MovieData.load(omeTiffPath);

% Set-up analysis infrastructure via command line interace
MD.addPackage(UTrackPackage(MD));
MD.getPackage(1).createDefaultProcess(1);

% Run the first process
MD.getPackage(1).getProcess(1).run();

% Setup the second process
MD.getPackage(1).createDefaultProcess(2);
MD.getPackage(1).getProcess(2).run();