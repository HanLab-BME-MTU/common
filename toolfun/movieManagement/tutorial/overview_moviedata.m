%% Create a MovieData object

%% MovieData initialization
% Initialize an MD from A the MAT file
dvPath = '~/Desktop/2014Mar20/110609_RhoWT_glycofect_001.dv';
MD = MovieData(dvPath);
fprintf(1, 'Object saved under: %s\n', MD.getFullPath());
fprintf(1, 'Output directory for analysis: %s\n', MD.outputDirectory_);

%% MovieData metadata

% Path location
disp('Raw data location');
fprintf(1, 'Channels path: %s\n', MD.getChannel(1).channelPath_);

% Retrieve movie dimensions
disp('Dimensions');
fprintf(1, '  Number of channels: %g\n', numel(MD.channels_));
fprintf(1, '  Number of timepoints: %g\n', MD.nFrames_);
fprintf(1, '  Image size: %gX%g\n', MD.imSize_);

% Retrieve raw metadata
disp('Metadata');
fprintf(1, '  Pixel size: %g nm\n', MD.pixelSize_);
fprintf(1, '  Numerical aperture: %g nm\n', MD.numAperture_);
for i = 1 : numel(MD.channels_),
    fprintf(1, '  Channel %g\n', i);
    fprintf(1, '    Emission wavelength: %g nm\n',...
        MD.getChannel(i).emissionWavelength_);
    fprintf(1, '    Excitation wavelength: %g nm\n',...
        MD.getChannel(i).excitationWavelength_);
    fprintf(1, '    Psf sigma: %g\n',...
        MD.getChannel(i).psfSigma_);
end

% Retrieve the raw data
disp('Pixel data');
for i = 1 : numel(MD.channels_)
    I = MD.getChannel(1).loadImage(1); % first timepoint of first channel
    fprintf(1, '  Plane 1 of channel %g of type %s and dimensions %gx%g\n',...
        i, class(I), size(I));
end

%% MovieData metadata

% Check
disp('Write Metadata');
fprintf(1, '  Initial time interval: %g s\n', MD.timeInterval_);

% Set time interval
MD.timeInterval_ = 1;
fprintf(1, '  Initial time interval: %g s\n', MD.timeInterval_);

% Reset same time interval
MD.timeInterval_ = 1;

% Setting different metadata fails
try
    MD.timeInterval_ = 2;
catch ME
    disp(ME.message)
end

% Setting invalid metadata fails
try
    MD.getChannel(1).exposureTime_ = - 2;
catch ME
    disp(ME.message)
end

%% Graphical User interface
% Launch the movie creation GUI
movieDataGUI();

% Launch the movie viewing GUI
movieViewer(MD);
