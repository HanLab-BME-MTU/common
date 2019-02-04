function [beadDensity]=calculateDensityFromDetection(MD)
%function []=calculateDensityFromDetection(MD) calculates bead density out
%of pointSourceDetection-based detections from the MD.
% Sangyoon Han September 2018

%% Get the beads number
% Load pointsource detection process
ipSD_Proc = MD.getProcessIndex('PointSourceDetectionProcess');
pSD_Proc = MD.getProcess(ipSD_Proc);

% Get the pointSource output from the process
pstruct = pSD_Proc.loadChannelOutput(1,1);

% Get the number
nBeads = length(pstruct.x);

% Get the channel
curChannel = MD.channels_(1);
curImage = curChannel.loadImage(1);

% Get the area
totalPixelArea = MD.imSize_(1) * MD.imSize_(2);
pixSize = MD.pixelSize_; %nm/pix
pixSizeMicron = pixSize/1000;
frameAreaMicron = totalPixelArea*pixSizeMicron^2;

% Get the density
beadDensity = nBeads/frameAreaMicron; % number/micron2

disp(['beadDensity= ' num2str(beadDensity) ' number/micron2.'])

% Save the density information to the folder where MD is stored
mkdir(MD.getPath, 'BeadDensity')
dataPath = [MD.getPath filesep 'BeadDensity'];
tableBeadDensity=table(beadDensity);
writetable(tableBeadDensity,[dataPath filesep 'beadDensity.csv'],'WriteRowNames',true)
