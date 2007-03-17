function gapInfo = findTrackGaps(trackedFeatureInfo)
%FINDTRACKGAPS finds the gaps in each track and gives back their location and length
%
%SYNOPSIS gapInfo = findTrackGaps(trackedFeatureInfo)
%
%INPUT  trackedFeatureInfo: Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%
%OUTPUT gapInfo           : An array with 3 columns. 1st column indicates
%                           track to which gap belongs, 2nd column 
%                           indicates time point where gap starts, 
%                           and 3rd column indicates gap length.
%
%Khuloud Jaqaman, February 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gapInfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--trackedFeatureInfo: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track information extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get track start times and end times
trackSEL = getTrackSEL(trackedFeatureInfo);

%make new matrix which contains only one column per time point
trackedFeatureInfo = trackedFeatureInfo(:,1:8:end);

%get number of tracks
numTracks = size(trackedFeatureInfo,1);

%alocate memory for output (assume that each track will have 10 gaps on average)
gapInfo = zeros(10*numTracks,3);

%assign value of index showing last stored position in gapInfo.
indxLast = 0;

%look for gaps in each track
for i=1:numTracks
    
    %get current track, its starting point and its ending point
    track0 = trackedFeatureInfo(i,:);
    start0 = trackSEL(i,1);
    end0 = trackSEL(i,2);
    
    %find all time points where track is NaN between its beginning and its end
    missing = find(isnan(track0))';
    missing = missing(missing>start0&missing<end0);

    %if there are gaps in track ...
    if ~isempty(missing)

        %find the time increment between one NaN and the next one
        missingDiff = diff(missing);

        %find all places where time increment is larger than 1 - this is the
        %start of a gap
        gapStart = missing([1;find(missingDiff~=1)+1]);

        %find places just before those with an increment larger than 1 - this
        %is the end of a gap
        gapEnd = missing([find(missingDiff~=1);end]);

        %calculate gap length
        gapLength = gapEnd - gapStart + 1;

        %get number of gaps in track
        numGaps = length(gapStart);

        %place the gaps of current track in gapInfo
        gapInfo(indxLast+1:indxLast+numGaps,:) = [i*ones(numGaps,1) gapStart gapLength];

        %update indxLast
        indxLast = indxLast + numGaps;

    end %(if ~isempty(missing))

end %(for i=1:numTracks)

%remove any unused entries in gapInfo
gapInfo = gapInfo(gapInfo(:,1)~=0,:);


%%%%% ~~ the end ~~ %%%%%

