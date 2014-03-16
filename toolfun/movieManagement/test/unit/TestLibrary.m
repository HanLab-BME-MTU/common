classdef TestLibrary < handle
    %TESTLIBRARY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        movie
        nRois
    end    
    
    methods
        
        function tearDown(self)
            cellfun(@delete, self.movie.processes_);
            cellfun(@delete, self.movie.packages_);
            delete(self.movie);
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
            rois(1, self.nRois) = MovieData();
            for i = 1 : self.nRois
                rois(i) = self.movie.addROI('','');
            end
        end
        
        function process = setUpProcess(self)
            process = MockProcess(self.movie);
            self.movie.addProcess(process);
        end
        
        function [package, process] = setUpPackage(self, loaded)
            package = MockPackage(self.movie);
            self.movie.addPackage(package);
            
            if nargin > 1 && loaded
                process = setUpProcess(self);
                package.setProcess(1, process);
            else
                process = [];
            end
        end
    end
end

