function stimulated = isMovieStimulated(movieData)

%
%Checks if there was a stimulation event recorded in the movieData. This
%could be a drug perfusion, serum withdrawl etc. which happens at some
%frame during the movie.
%
%Hunter Elliott 3/2009
%

%Check that the stimulation information is there and makes sense
if isfield(movieData,'stimulation') && isfield(movieData.stimulation,'iFrame') ...
        && ~isempty(movieData.stimulation.iFrame) && movieData.stimulation.iFrame > 0 ...
        && movieData.stimulation.iFrame == round(movieData.stimulation.iFrame)
        
    stimulated = true;
else
    stimulated = false;
end