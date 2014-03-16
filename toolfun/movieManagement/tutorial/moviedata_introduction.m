%% MOVIEDATA introduction

% Initialize an empty MovieData object
MD = MovieData();

% Check class properties
disp(class(MD)); % Display the class
disp(isa(MD, 'MovieData')); 
disp(isa(MD, 'MovieObject'));
disp(isa(MD, 'handle'));

% Display the object properties and methods
properties(MD);
methods(MD);

%% Documentation

% Use MATLAB built-in help
help MovieData
help MovieData.camBitdepth_
help MovieData.load

% Use MATLAB built-in doc
doc MovieData
doc MovieData.camBitDepth


%% Metadata

