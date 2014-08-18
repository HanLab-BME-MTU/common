function movieInfo = correctIntensity(movieInfo,correctionImage)
%CORRECTINTENSITY corrects object intensities due to uneven illumination
%
%SYNOPSIS movieInfo = correctMovieInfoInt(movieInfo,correctionImage)
%
%INPUT  
%       movieInfo    : Array of size equal to the number of frames
%                      in movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%             .yCoord      : Image coordinate system y-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D. Default: zeros.
%             .zCoord      : Image coordinate system z-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D or 2D. Default: zeros.
%             .amp         : Amplitudes of PSFs fitting detected features. 
%                            1st column for values and 2nd column 
%                            for standard deviations.
%       correctionImage: Image to use for correction.
%
%OUTPUT movieInfo    : Same as input, but with added field "ampUncorr,"
%                      where the original uncorrected intensities are moved
%                      to the field "ampUncorr" and the corrected amplitudes
%                      are saved in the old field "amp."
%
%Khuloud Jaqaman, August 2014

%% Input

%find number of frames in movie
numFrames = length(movieInfo);

%find image size
correctionImage = double(correctionImage);
imSize = size(correctionImage);

%% Correction

%normalize correction image so that average is 1
meanInt = mean(correctionImage(:));
correctionImage = correctionImage / meanInt;

%go over all frames
for iFrame = 1 : numFrames
   
    %move uncorrected amplitudes to new field
    movieInfo(iFrame).ampUncorr = movieInfo(iFrame).amp;
    
    if ~isempty(movieInfo(iFrame).xCoord)
        
        %object information
        xCoord = round(movieInfo(iFrame).xCoord(:,1));
        yCoord = round(movieInfo(iFrame).yCoord(:,1));
        amp = movieInfo(iFrame).amp(:,1);
        
        %linear index
        linIdx = sub2ind(imSize,yCoord,xCoord);
        
        %correction image value at linIdx
        correctionVal = correctionImage(linIdx);
        
        %correct amplitude
        amp = amp ./ correctionVal;
        movieInfo(iFrame).amp(:,1) = amp;
        
    end
    
end

%%%%% ~~ the end ~~ %%%%%
