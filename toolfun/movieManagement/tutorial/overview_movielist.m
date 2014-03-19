%% Create a MovieList object

clear
clc
close all

%% Initialization/constructor
% Define original paths
rootPath = '~/Desktop/2014Mar20/MovieList';
movieNames = {'P4A', 'P14', 'P15'};

% Load individual movies
MD(numel(movieNames), 1)=  MovieData;
for i = 1 : numel(movieNames)
    MD(i) = MovieData.load(fullfile(rootPath, movieNames{i}, 'movieData.mat'));
end

% Create movie list
ML = MovieList(MD, rootPath);

%% Manipulation via CLI

% Set path properties
ML.setPath(rootPath);
ML.setFilename('movieList.mat');

% Save movie
ML.save();

%% Movie access
% Retrieve individual movies
disp('Movies')
for i = 1 : numel(ML.movieDataFile_)
    fprintf(1, '  Movie %g: %s\n', i, ML.getMovie(i).getFullPath());
end

%
disp('Output')
fprintf(1, '  Analysis saved under: %s\n', ML.outputDirectory_);

%% Graphical interface

% Launch viewing interface
movieViewer(ML);