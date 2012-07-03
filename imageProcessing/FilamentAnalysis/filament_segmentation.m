function movieData = filament_segmentation(movieData, varargin)

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('funParams',[], @isstruct);
ip.parse(movieData,varargin{:});
funParams=ip.Results.funParams;

