%makeMovieFFMPEG(stack, varargin) generates a movie from the input stack using FFMPEG and libx264.
% FFMPEG with libx264 must be installed, and this function only works on a Unix system.
%
% Inputs: 
%        stack : 3D stack of movie frames
%
% Options:
%
%    framerate : the frame rate of the output movie. Default: 15
%
% Parameters:
%   'Destpath' : destination directory for the movie. Default: current directory (pwd)
%    'Quality' : Quality setting for ffmpeg. Default: 22

% Francois Aguet, 10/16/2013

function makeMovieFFMPEG(stack, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('stack');
ip.addOptional('framerate', 15, @isposint);
ip.addParamValue('DestPath', []); % for movie
ip.addParamValue('FileName', 'movie.mp4', @ischar);
ip.addParamValue('Quality', 22, @isposint);
ip.parse(stack, varargin{:});

[ny,nx,nf] = size(stack);
fmt = ['%0' num2str(ceil(log10(nf))) 'd'];

if isunix && ~system('which ffmpeg >/dev/null 2>&1')
    fprintf('Generating movie ... ');
    frameDest = ['.frames_tmp_mmffmpeg' filesep];
    [~,~] = mkdir(frameDest);
    dRange = double([min(stack(:)) max(stack(:))]);
    for fi = 1:nf
        imwrite(uint8(scaleContrast(stack(:,:,fi), dRange)), [frameDest 'frame_' num2str(fi, fmt), '.png']);
    end
    
    fr = num2str(ip.Results.framerate);
    
    cmd = ['ffmpeg -quiet -y -r ' fr ' -i ' frameDest 'frame_' fmt '.png' ' -vf "scale=' num2str(2*floor(nx/2)) ':' num2str(2*floor(ny/2))...
        '" -c:v libx264 -crf ' num2str(ip.Results.Quality) ' -pix_fmt yuv420p ' ip.Results.DestPath ip.Results.FileName];
    system(cmd);
    rmdir(frameDest, 's');
    fprintf(' done.\n');
else
    fprintf('A unix system with ffmpeg installed is required to generate movies using this function.\n');
end
