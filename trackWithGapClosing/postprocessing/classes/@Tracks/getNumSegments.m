function [ N ] = getNumSegments( obj , idx)
%numSegments gets the number of segments within an array of compound tracks
if(nargin < 2)
    tracks = obj;
else
    tracks = obj(idx);
end

    N = [tracks.numSegments];

end

