function movie2aviNADA(mov, filename, varargin)

% movie2aviNADA(mov, filename, ...)
%
% movie2avi, Not A DumbAss version
%
% since movie2avi barfs if any frames of mov are empty, I had to write this
% function which deals with that first and then calls movie2avi on the
% modified mov; currently black frames are inserted
%
% Cyrus A Wilson    July 2006

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

numFrames = max(size(mov));

%% first go through and find first non-empty frame

for m = 1:numFrames
  if ~isempty(mov(m).cdata)
    defcdata = mov(m).cdata;
    defcolormap = mov(m).colormap;
    break;
  end
end

if ~exist('defcdata', 'var')
  % no non-empty frames; nothing to do so just return
  return;
end

defcdata(:) = 0; %set to black (could cause problems for indexed images)

%% now fill in empty frames

for m = 1:numFrames
  if isempty(mov(m).cdata)
    mov(m).cdata = defcdata;
    mov(m).colormap = defcolormap;
  end
end

%% now call movie2avi on the modified mov

movie2avi(mov, filename, varargin{:});
