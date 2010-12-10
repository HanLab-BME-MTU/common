function [movieData,intStats] = calculateMovieIntensityVsTime(movieData,p)
%CALCULATEMOVIEINTENSITYVSTIME calculates the intensity in each frame of the input movie
%
% movieData = calculateMovieIntensityVsTime(movieData)
% [movieData,intensityStatistics] = calculateMovieIntensityVsTime(movieData)
% ... = calculateMovieIntensityVsTime(movieData,paramIn)
%
%   This function calculates various statistics (mean, STD, median etc.)
%   regarding the pixel intensity values in each frame of the selected
%   channel(s) of the input movie. This can be either the masked intensity
%   or the intensity of the entire image. The results are saved to the
%   movie's output directory and returned as output. A figure is also
%   created and saved to this directory.
% 
% Input:
% 
%   movieData - A MovieData object describing the movie to be processed, as
%   created by setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       results to. Optional. If not specified, the results will be saved
%       to the movieData's output directory in a sub-folder called
%       "intensityVsTime"
%
%       ('ChannelIndex' -> Positive integer scalar or vector) Optional. The
%       integer index of the channel(s) to calculate statistics for. If not
%       input, all channels will be used.
% 
%       ('MaskedIntensity'-> True/False) If true, the statistics will be
%       calculated only for masked pixels. This requires that the movie has
%       already been segmented. If false, all pixels are used. Optional.
%       Default is false.
%
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes, figures etc.)
%       Optional. Default is false.
%
% Output:
% 
%   movieData - The updated movieData object with the calculations logged
%   in the processes_ array.
% 
%   intensityStatistics - An Mx1 array of structures each of which contains
%   the following fields:
%       
% 
% 
% 
% 
% Hunter Elliott 
% 11/2010
%


%TEMP - STILL UNDER CONSTRUCTION

fName = 'masked_intensity_vs_time';
dName = 'intensityVsTime';

p.BatchMode = false;

p.SegProcessIndex = movieData.getProcessIndex('MaskRefinementProcess',1,1);

p.ProcessIndex = movieData.getProcessIndex('BleedthroughCorrectionProcess',1,1);

p.ProcessIndex = [p.ProcessIndex movieData.getProcessIndex('BackgroundSubtractionProcess',1,1)];

if isempty(p.SegProcessIndex)
    error('No valid segmentation process found!')
end

if isempty(p.MaskChannelIndex)
    p.MaskChannelIndex = p.ChannelIndex;
end


nChan = numel(p.ChannelIndex);
if nChan == 0
    error('You must specify at least one channel!');
end

for j = 1:nChan
    
    for k = 1:numel(p.ProcessIndex)
        if movieData.processes_{p.ProcessIndex(k)}.checkChannelOutput(p.ChannelIndex(j));
            imNames(j) = movieData.processes_{p.ProcessIndex(k)}.getOutImageFileNames(p.ChannelIndex(j));            
            imDirs{j} = movieData.processes_{p.ProcessIndex(k)}.outImagePaths_{p.ChannelIndex(j)};
            break
        end
        if k == numel(p.ProcessIndex)
            error(['Could not find valid images for channel ' num2str(p.ChannelIndex(j))])
        end
        
    end
end

maskDirs = movieData.processes_{p.SegProcessIndex}.outMaskPaths_(p.MaskChannelIndex);
maskNames = movieData.processes_{p.SegProcessIndex}.getOutMaskFileNames(p.MaskChannelIndex);


p.OutputDirectory = [movieData.outputDirectory_ filesep dName];
mkClrDir(p.OutputDirectory);


nImages = movieData.nFrames_;


avgIntensity = nan(nChan,nImages);
stdIntensity = nan(nChan,nImages);
nPix = nan(nChan,nImages);
for iChan = 1:nChan
    
    for iImage = 1:nImages
        
        
        currImage = imread([imDirs{iChan} filesep imNames{iChan}{iImage}]);
        
        currMask = imread([maskDirs{iChan} filesep maskNames{iChan}{iImage}]);
        
        avgIntensity(iChan,iImage) = mean(currImage(currMask(:)));
        stdIntensity(iChan,iImage) = std(double(currImage(currMask(:))));
        nPix(iChan,iImage) = nnz(currMask);
        
    end
    
end

save([p.OutputDirectory filesep fName '.mat'],'avgIntensity',...
                                              'stdIntensity',...
                                              'nPix');

if ~p.BatchMode
    chanCols = lines(nChan);
    figure    
    hold on
    for j = 1:nChan
        plot(avgIntensity(j,:),'color',chanCols(j,:));
        leg{j} = ['Channel ' num2str(p.ChannelIndex(j)) ];
        %plotTransparent(1:nImages,avgIntensity(j,:),stdIntensity(j,:),chanCols(j,:),.2,0);
    end
    legend(leg{:});
end
    



