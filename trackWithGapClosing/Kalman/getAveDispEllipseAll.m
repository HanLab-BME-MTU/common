function [dispDrift,dispBrown,longVec,shortVec] = getAveDispEllipseAll(xVel,...
    yVel,brownStd,trackType,undetBrownStd,timeWindow,brownStdMult,linStdMult,...
    timeReachConf,minSearchRadius,useLocalDensity,closestDistScale,...
    maxStdMult,nnDistTracks)
%GETAVEDISPELLIPSE determines the search ellipse and expected displacement along x and y of a particle undergoing 2D diffusion with drift
%
%SYNOPSIS [dispDrift,dispBrown,longVec,shortVec] = getAveDispEllipseAll(xVel,...
%    yVel,brownStd,trackType,undetBrownStd,timeWindow,brownStdMult,linStdMult,...
%    timeReachConf,minSearchRadius,useLocalDensity,closestDistScale,maxStdMult)
%
%INPUT  xVel           : Velocity in x-direction.
%       yVel           : Velocity in y-direction.
%       brownStd       : Standard deviation of Brownian motion steps.
%       trackType      : Type of track. 1 for directed, 0 for Brownian, NaN for unetermined.
%       undetBrownStd  : Standard deviation of Brownian motion steps to be used
%                        for undetermined tracks.
%       timeWindow     : Maximum gap size.
%       brownStdMult   : Multiplication factor to go from average Brownian
%                        displacement to search radius.
%       linStdMult     : Multiplication factor to go from average linear
%                        displacement to search radius.
%       timeReachConf  : Time gap for reaching confinement.
%       minSearchRadius: Minimum allowed search radius.
%       useLocalDensity: 1 if local density of features is used to expand 
%                        their search radius if possible, 0 otherwise.
%       closestDistScale:Scaling factor of nearest neighbor distance.
%       maxStdMult     : Maximum value of factor multiplying std to get
%                        search radius.
%       nnDistTracks   : Vector indicating the nearest neighbor distance
%                        of each track, i.e. the closest distance any of
%                        the features making it comes to any other
%                        feature in the feature's frame.
%
%OUTPUT dispDrift : Column vector of expected displacement along x and y due to drift.
%       dispBrown : Expected displacement along x (= along y) due to
%                   Brownian motion.
%       longVec   : Column vector defining long radius of search ellipse.
%       shortVec  : Column vector defining short radius of search ellipse.
%       errFlag   : 0 if function executes normally, 1 otherwise
%
%REMARKS Drift is assumed to look more like 1D diffusion, i.e. the particle
%goes back and forth along a line
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dispDrift = [];
dispBrown = [];
longVec = [];
shortVec = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('getAveDispEllipseAll')
    disp('--getAveDispEllipseAll: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine expected displacement and search ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine number of tracks
numTracks = length(xVel);

%reserve memory for output
dispDrift = zeros(2,timeWindow,numTracks);
dispBrown = zeros(timeWindow,numTracks);
longVec = zeros(2,timeWindow,numTracks);
shortVec = zeros(2,timeWindow,numTracks);

%define square root of two to avoid calculating it many times
sqrtTwo = sqrt(2);

%put time scaling of linear motion in a vector
sqrtTimeGap = sqrt(1:timeWindow);

%put time scaling of Brownian motion in a vector
timeScalingBrown = [sqrt(1:timeReachConf) sqrt(timeReachConf) * ...
    (2:timeWindow-timeReachConf+1).^0.1];

for iTrack = 1 : numTracks

    switch trackType(iTrack)

        case 1
            
            %get velocity, its magnitude and direction of motion
            velDrift = [xVel(iTrack) yVel(iTrack)]';
            velMag = sqrt(velDrift' * velDrift);
            directionMotion = velDrift / velMag;

            %calculate the expected displacement due to drift for all time
            %gaps
            dispDrift1 = repmat(velDrift,1,timeWindow) .* repmat(sqrtTimeGap,2,1);

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;

            %determine the long vector of the search ellipse for all time
            %gaps
            longVec1 = repmat(linStdMult',2,1) .* dispDrift1 + ...
                repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);
            
            %determine the short vector
            shortVec1 = repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %output the absolute value of dispDrift
            dispDrift1 = abs(dispDrift1);

            %make sure that long vectors are longer than minimum allowed
            vecMag = sqrt((diag(longVec1' * longVec1))');  %magnitude
            vecDir = longVec1 ./ repmat(vecMag,2,1); %direction
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1 = repmat(vecMag,2,1) .* vecDir; %new long vector

            %make sure that short vectors are longer than minimum allowed
            vecMag = sqrt((diag(shortVec1' * shortVec1))');  %magnitude
            vecDir = shortVec1 ./ repmat(vecMag,2,1); %direction
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVec1 = repmat(vecMag,2,1) .* vecDir; %new short vector

            %save values for this time gap
            dispDrift(:,:,iTrack) = dispDrift1;
            dispBrown(:,iTrack) = dispBrown1;
            longVec(:,:,iTrack) = longVec1;
            shortVec(:,:,iTrack) = shortVec1;

        case 0
            
            %take direction of motion along x
            directionMotion = [1 0]';

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;

            %determine the long vector of the search ellipse for all time
            %gaps
            longVec1 = repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vector
            shortVec1 = repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %make sure that long vectors are longer than minimum allowed
            vecMag = sqrt((diag(longVec1' * longVec1))');  %magnitude
            vecDir = longVec1 ./ repmat(vecMag,2,1); %direction
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1 = repmat(vecMag,2,1) .* vecDir; %new long vector

            %make sure that short vectors are longer than minimum allowed
            vecMag = sqrt((diag(shortVec1' * shortVec1))');  %magnitude
            vecDir = shortVec1 ./ repmat(vecMag,2,1); %direction
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVec1 = repmat(vecMag,2,1) .* vecDir; %new short vector

            %save values for this time gap
            dispBrown(:,iTrack) = dispBrown1;
            longVec(:,:,iTrack) = longVec1;
            shortVec(:,:,iTrack) = shortVec1;

        otherwise

            %take direction of motion along x
            directionMotion = [1 0]';

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = undetBrownStd * timeScalingBrown;

            %determine the long vector of the search ellipse for all time
            %gaps
            longVec1 = repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vector
            shortVec1 = repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %make sure that long vectors are longer than minimum allowed
            vecMag = sqrt((diag(longVec1' * longVec1))');  %magnitude
            vecDir = longVec1 ./ repmat(vecMag,2,1); %direction
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1 = repmat(vecMag,2,1) .* vecDir; %new long vector

            %make sure that short vectors are longer than minimum allowed
            vecMag = sqrt((diag(shortVec1' * shortVec1))');  %magnitude
            vecDir = shortVec1 ./ repmat(vecMag,2,1); %direction
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVec1 = repmat(vecMag,2,1) .* vecDir; %new short vector

            %save values for this time gap
            dispBrown(:,iTrack) = dispBrown1;
            longVec(:,:,iTrack) = longVec1;
            shortVec(:,:,iTrack) = shortVec1;

    end %(switch trackType)

end %(for iTrack = 1 : numTracks)

%%%%% ~~ the end ~~ %%%%%

