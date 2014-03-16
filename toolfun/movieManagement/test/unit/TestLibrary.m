classdef TestLibrary < handle
    %TESTLIBRARY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        movie
        nRois
        rois = MovieData.empty(0,1);
        packages = {};
        processes = {};
    end    
    
    methods
        
        function tearDown(self)
            delete(self.movie);
            cellfun(@delete, self.processes);
            cellfun(@delete, self.packages);
        end
        
        function setUpMovieData(self)
            self.movie = MovieData();
        end
        
        
        function setUpMovieList(self)
            self.movie = MovieList();
        end
        
        function rois = setUpRois(self, varargin)
            if nargin > 1
                self.nRois = varargin{1};
            else
                self.nRois  = 5;
            end
            for i = 1 : self.nRois
                self.rois(end+1) = self.movie.addROI('','');
            end
            rois = self.rois(end - self.nRois + 1 : end);
        end
        
        function process = setUpProcess(self)
            process = MockProcess(self.movie);
            self.movie.addProcess(process);
            self.processes{end+1} = process;
        end
        
        function [package, process] = setUpPackage(self, loaded)
            package = MockPackage(self.movie);
            self.movie.addPackage(package);
            self.packages{end+1} = package;
            
            if nargin > 1 && loaded
                process = setUpProcess(self);
                package.setProcess(1, process);
            else
                process = [];
            end
        end
    end
end

