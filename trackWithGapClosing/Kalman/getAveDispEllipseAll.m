function [dispDrift,dispBrown,longVecS,longVecE,shortVecS,shortVecE] = ...
    getAveDispEllipseAll(xVel,yVel,brownStd,trackType,undetBrownStd,...
    timeWindow,brownStdMult,linStdMult,timeReachConfB,timeReachConfL,...
    minSearchRadius,maxSearchRadius,useLocalDensity,closestDistScale,...
    maxStdMult,nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime)
%GETAVEDISPELLIPSE determines the search ellipse and expected displacement along x and y of a particle undergoing 2D diffusion with drift
%
%SYNOPSIS [dispDrift,dispBrown,longVecS,longVecE,shortVecS,shortVecE] = ...
%    getAveDispEllipseAll(xVel,yVel,brownStd,trackType,undetBrownStd,...
%    timeWindow,brownStdMult,linStdMult,timeReachConfB,timeReachConfL,...
%    minSearchRadius,maxSearchRadius,useLocalDensity,closestDistScale,...
%    maxStdMult,nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime)
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
%       timeReachConfB : Time gap for Brownian motion to reach confinement.
%       timeReachConfL : Time gap for linear motion to reach confinement.
%       minSearchRadius: Minimum allowed search radius.
%       maxSearchRadius: Maximum allowed search radius for linking between
%                        two consecutive frames. It will be expanded for
%                        different gap lengths based on the time scaling of
%                        Brownian motion.
%       useLocalDensity: 1 if local density of features is used to expand 
%                        their search radius if possible, 0 otherwise.
%       closestDistScale:Scaling factor of nearest neighbor distance.
%       maxStdMult     : Maximum value of factor multiplying std to get
%                        search radius.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       nnWindow       : Time window to be used in estimating the
%                        nearest neighbor distance of a track at its start
%                        and end.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%
%OUTPUT dispDrift : Vector of expected displacement along x and y due to drift.
%       dispBrown : Expected displacement along x (= along y) due to
%                   Brownian motion.
%       longVecS  : Vector defining long radius of search ellipse at the
%                   starts of tracks.
%       longVecE  : Vector defining long radius of search ellipse at the
%                   ends of tracks.
%       shortVecS : Vector defining short radius of search ellipse at the
%                   starts of tracks.
%       shortvecE : Vector defining short radius of search ellipse at the
%                   ends of tracks.
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
longVecS  = [];
longVecE  = [];
shortVecS = [];
shortVecE = [];

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
longVecS  = zeros(2,timeWindow,numTracks);
longVecE  = zeros(2,timeWindow,numTracks);
shortVecS = zeros(2,timeWindow,numTracks);
shortVecE = zeros(2,timeWindow,numTracks);

%define square root of two to avoid calculating it many times
sqrtTwo = sqrt(2);

%put time scaling of linear motion in a vector
timeScalingLin = [sqrt(1:timeReachConfL) sqrt(timeReachConfL) * ...
    (2:timeWindow-timeReachConfL+1).^0.1];

%put time scaling of Brownian motion in a vector
timeScalingBrown = [sqrt(1:timeReachConfB) sqrt(timeReachConfB) * ...
    (2:timeWindow-timeReachConfB+1).^0.1];

%scale maxSearchRadius like Brownian motion (it's only imposed on the
%Brownian aspect of tracks)
maxSearchRadius = maxSearchRadius * timeScalingBrown;

%determine the nearest neighbor distances of tracks at their starts and ends
windowLimS = min([trackStartTime+nnWindow trackEndTime],[],2);
windowLimE = max([trackEndTime-nnWindow trackStartTime],[],2);
nnDistTracksS = zeros(numTracks,1);
nnDistTracksE = zeros(numTracks,1);
for iTrack = 1 : numTracks
    nnDistTracksS(iTrack) = min(nnDistLinkedFeat(iTrack,...
        trackStartTime(iTrack):windowLimS(iTrack)));
    nnDistTracksE(iTrack) = min(nnDistLinkedFeat(iTrack,...
        windowLimE(iTrack):trackEndTime(iTrack)));
end

for iTrack = 1 : numTracks

    switch trackType(iTrack)

        case 1
            
            %get velocity, its magnitude and direction of motion
            velDrift = [xVel(iTrack) yVel(iTrack)]';
            velMag = sqrt(velDrift' * velDrift);
            directionMotion = velDrift / velMag;

            %calculate the expected displacement due to drift for all time
            %gaps
            dispDrift1 = repmat(velDrift,1,timeWindow) .* repmat(timeScalingLin,2,1);

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;
            
            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance at its start
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track start
                %if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track end
                %if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the long vectors of the search ellipses for all time
            %gaps
            longVec1 = repmat(linStdMult',2,1) .* dispDrift1 + ...
                repmat((brownStdMult' .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);
            longVecMag = sqrt((diag(longVec1' * longVec1))');  %magnitude
            longVecDir = longVec1 ./ repmat(longVecMag,2,1); %direction
            
            %determine the short vectors at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);
            shortVecSMag = sqrt((diag(shortVecS1' * shortVecS1))');  %magnitude
            shortVecSDir = shortVecS1 ./ repmat(shortVecSMag,2,1); %direction

            %determine the short vectors at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);
            shortVecEMag = sqrt((diag(shortVecE1' * shortVecE1))');  %magnitude
            shortVecEDir = shortVecE1 ./ repmat(shortVecEMag,2,1); %direction

            %output the absolute value of dispDrift
            dispDrift1 = abs(dispDrift1);

            %make sure that long vectors are longer than minimum allowed
            longVecMag = max([longVecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1 = repmat(longVecMag,2,1) .* longVecDir; %new long vector

            %make sure that short vectors at track starts are longer than 
            %minimum allowed and shorter than maximum allowed
            shortVecSMag = max([shortVecSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVecSMag = min([shortVecSMag;maxSearchRadius]); %compare to maximum
            shortVecS1 = repmat(shortVecSMag,2,1) .* shortVecSDir; %new short vector

            %make sure that short vectors at track ends are longer than 
            %minimum allowed and shorter than maximum allowed
            shortVecEMag = max([shortVecEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVecEMag = min([shortVecEMag;maxSearchRadius]); %compare to maximum
            shortVecE1 = repmat(shortVecEMag,2,1) .* shortVecEDir; %new short vector

            %save values for this time gap
            dispDrift(:,:,iTrack) = dispDrift1;
            dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVec1;
            longVecE(:,:,iTrack) = longVec1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;

        case 0
            
            %take direction of motion to be along x
            directionMotion = [1 0]';

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end
            
            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track end
                %if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the long vectors of the search ellipses at track
            %starts for all time gaps
            longVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vectors of the search ellipses at track
            %ends for all time gaps
            longVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vectors at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %determine the short vectors at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,2,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,2,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecS1 = repmat(vecMag,2,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMag,2,1) .* shortVecDir; %new short vector

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))');  %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,2,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,2,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecE1 = repmat(vecMag,2,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMag,2,1) .* shortVecDir; %new short vector

            %save values for this time gap
            dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;

        otherwise

            %take direction of motion to be along x
            directionMotion = [1 0]';

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = undetBrownStd * timeScalingBrown;

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end
            
            %determine the long vector of the search ellipse at track
            %starts for all time gaps
            longVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vector of the search ellipse at track
            %ends for all time gaps
            longVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vector at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %determine the short vector at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtTwo),2,1) .* ...
                repmat([directionMotion(2) -directionMotion(1)]',1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,2,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,2,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecS1 = repmat(vecMag,2,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMag,2,1) .* shortVecDir; %new short vector

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))'); %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,2,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,2,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecE1 = repmat(vecMag,2,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMag,2,1) .* shortVecDir; %new short vector

            %save values for this time gap
            dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;

    end %(switch trackType)

end %(for iTrack = 1 : numTracks)

%%%%% ~~ the end ~~ %%%%%

