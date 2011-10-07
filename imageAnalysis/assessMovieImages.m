function assessMovieImages(movieData,varargin)
%ASSESSMOVIEIMAGEQUALITY performs basic image quality assessment on the images in the input movie 
% 
% assessMovieImageQuality(movieData)
% assessMovieImageQuality(movieData,'OptionName',optionValue,...)
% 
% This function performs some calculations to assess the quality of the
% fluorescence images in the input movie. Various figures showing measures
% of the image quality will be displayed, and the statistics will be saved
% to the specified output directory.
% Note that some properties require that masks have already been created
% for the input movie. These are indicated by an asterix below.
%
%   Calculated properties include (* indicates masks required):
%       -image intensity statistics over time (avg, min, max, etc.)
%       -photobleaching rate
%       -masked/unmasked intensity statistics*
%       -foreground/backround intensity statistics*
%       -signal:background ratio*
%       -saturation amount 
%        
% 
% Input:
%
%   movieData - A MovieData object describing the movie to be processed, as
%   created by setupMovieDataGUI.m
%
%   'OptionName',optionValue - An option name character array followed by
%   the value for that option. Options are described below.
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       statistics and figures to.
%       If not input, the output will be saved to the same directory as the
%       movieData, in a sub-directory called "image_assessment"
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       Optional. The integer index of the channel(s) to assess. If not
%       input, all channels will be assessed.
%
%       ('SegProcessIndex' -> Positive integer scalar) This specifies the
%       index of the segmentation process to use masks from for calculating
%       foreground/background statistics etc. This is only necessary if more
%       than one segmentatin process exists. If not specified, and more
%       than one segmentation process exists, the user will be asked,
%       unless batch mode is enabled, in which case an error will be
%       generated.
%
%       ('BackSegProcessIndex' -> Positive integer scalar) This specifies
%       the index of the background masking process to use masks from for
%       calculating background image statistics. This is only
%       necessary if more than one segmentatin process exists. If not
%       specified, and more than one segmentation process exists, the user
%       will be asked, unless batch mode is enabled, in which case an error
%       will be generated.
%
%       ('MaskChannelIndex' -> Positive integer scalar or vector)
%       This parameter specifies the channels to use masks from when
%       calculating foreground/background statistics etc. This allows
%       masks generated for one channel to be used in performing background
%       subtraction on other channels, which may or may not themselves have
%       masks. This vector or scalar must be the same size as ChannelIndex.
%       Optional. If not input, each channel will use it's own masks.
% 
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
% 
% Output:
%
%   The statistics and figures will be saved to disk in the specified
%   output directory.
%  
%   If not in batch mode, figures will be displayed with the image
%   assessment statistics.
%
%
% Hunter Elliott
% 5/2011
%

%% ----------------------- Parameters ---------------------- %%

outName = 'image_assessment';%Sub-dir of movie outputDirectory to save results in if not specified
%Print options for saving figures to disk
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));%Double-exponential function for fitting bleaching curves
bInit = [1 0 1 0]; %Initial guess for fit parameters.
fitOptions = statset('Robust','on','MaxIter',1e3,'Display','off','FunValCheck','off');%Settings for fitting with nlinfit

%% ------------------------- Input -------------------------- %%

ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('movieData',@(x)(isa(x,'MovieData')));
ip.addParamValue('OutputDirectory',[],@ischar);
ip.addParamValue('ChannelIndex',[],@(x)(all(isposint(x))));
ip.addParamValue('SegProcessIndex',[],@(x)(isposint(x)));
ip.addParamValue('BackSegProcessIndex',[],@(x)(isposint(x)));
ip.addParamValue('MaskChannelIndex',[],@(x)(all(isposint(x))));
ip.addParamValue('BatchMode',false,(@(x)(numel(x)==1)));
ip.parse(movieData,varargin{:});

p = ip.Results;

if isempty(p.OutputDirectory)
    p.OutputDirectory = [movieData.outputDirectory_ filesep outName];
    mkClrDir(p.OutputDirectory);
elseif ~exist(p.OutputDirectory,'dir')
    error('The specified output directory does not exist!')
end

if isempty(p.ChannelIndex)
    if ~p.BatchMode
        p.ChannelIndex = selectMovieChannels(movieData,1,'Select the channels to assess:');
    else
        p.ChannelIndex = 1:numel(movieData.channels_);
    end
elseif any(p.ChannelIndex > numel(movieData.channels_))
    error('Invalid channel indices specified! Check ChannelIndex option!')
end

if isempty(p.SegProcessIndex)
    if ~p.BatchMode
        p.SegProcessIndex = movieData.getProcessIndex('MaskProcess',1,1);
    else
        p.SegProcessIndex = movieData.getProcessIndex('MaskProcess',Inf,0);
        if numel(p.SegProcessIndex) > 1
            error('If batch mode is not enabled and the SegProcessIndex is not specified, the movie cannot have more than one MaskProcess!')
        end
    end
end

if isempty(p.MaskChannelIndex)
    p.MaskChannelIndex = p.ChannelIndex;
elseif any(p.MaskChannelIndex) > numel(movieData.channels_) || ...
    numel(p.MaskChannelIndex) ~= numel(p.ChannelIndex)
    error('Invalid MaskChannelIndex : you must specify one valid mask channel index for every channel which is assessed!')
end

if isempty(p.SegProcessIndex) || ~isa(movieData.processes_{p.SegProcessIndex},'MaskProcess') ...
        || ~all(movieData.processes_{p.SegProcessIndex}.checkChannelOutput(p.MaskChannelIndex))
    warning('ImageAssessment:noMasks','A valid segmentation process for the specified channels was not present and specified by the user: some statistics will not be calculated!')
    hasMasks = false;
else
    hasMasks = true;
end

if isempty(p.BackSegProcessIndex)
    if ~p.BatchMode
        p.BackSegProcessIndex = movieData.getProcessIndex('BackgroundMasksProcess',1,1);
    else
        p.BackSegProcessIndex = movieData.getProcessIndex('BackgroundMasksProcess',Inf,0);
        if numel(p.BackSegProcessIndex) > 1
            error('If batch mode is not enabled and the SegProcessIndex is not specified, the movie cannot have more than one MaskProcess!')
        end
    end
end

if isempty(p.BackSegProcessIndex) || ~isa(movieData.processes_{p.BackSegProcessIndex},'BackgroundMasksProcess') ...
        || ~all(movieData.processes_{p.BackSegProcessIndex}.checkChannelOutput(p.MaskChannelIndex))
    warning('ImageAssessment:noMasks','A valid background segmentation process for the specified channels was not present and specified by the user: some statistics will not be calculated!')
    hasBackMasks = false;
else
    hasBackMasks = true;
end

if hasBackMasks && ~hasMasks
    error('Background masks found but no foreground masks available/selected! You must have masks if you have background masks!')
end

%% ------------------------ Init ----------------------- %%

imDirs = movieData.getChannelPaths(p.ChannelIndex);
imNames = movieData.getImageFileNames(p.ChannelIndex);


%If available, get mask directories and file names
if hasMasks
    maskDirs = movieData.processes_{p.SegProcessIndex}.outFilePaths_(p.MaskChannelIndex);
    maskNames = movieData.processes_{p.SegProcessIndex}.getOutMaskFileNames(p.MaskChannelIndex);
end

if hasBackMasks
    backMaskDirs = movieData.processes_{p.BackSegProcessIndex}.outFilePaths_(p.MaskChannelIndex);
    backMaskNames = movieData.processes_{p.BackSegProcessIndex}.getOutMaskFileNames(p.MaskChannelIndex);
end

nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
nImTot = nFrames * nChan;

if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, assessing channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end        
    
%Initialize structure for storing stats
statFields = {'AverageIntensity',...
              'MinimumIntensity',...
              'MaximumIntensity',...
              'MedianIntensity',...
              'STDIntensity',...
              'BleachingStats'};          
if hasMasks              
    statFields = [statFields ...
               {'AverageMaskedIntensity',...
                'MinimumMaskedIntensity',...
                'MaximumMaskedIntensity',...
                'MedianMaskedIntensity',...
                'STDMaskedIntensity',...
                'AverageUnmaskedIntensity',...
                'MinimumUnmaskedIntensity',...
                'MaximumUnmaskedIntensity',...
                'MedianUnmaskedIntensity',...
                'STDUnmaskedIntensity',...      
                'MaskedBleachingStats',...
                'RatioMaskedToUnmaskedIntensity'}];
end
if hasBackMasks
    
    statFields = [statFields ...
                {'AverageBackgroundIntensity',...
                 'MinimumBackgroundIntensity',...
                 'MaximumBackgroundIntensity',...
                 'MedianBackgroundIntensity',...
                 'STDBackgroundIntensity',...
                 'RatioMaskedToBackgroundIntensity'}];
end

nFields = numel(statFields);
imStats(1:nChan) = struct;
for j = 1:nChan
    for k = 1:nFields
        imStats(j).(statFields{k}) = nan(1,nFrames);
    end
end

%Since we overwrite the BleachStats field, disable the warning
warning('off','MATLAB:warn_r14_stucture_assignment')

%Options for figure creation/display
if p.BatchMode
    figOpts = {'Visible','off'};
else
    figOpts = {};
end
    
%Time data, if available
if ~isempty(movieData.timeInterval_)
    tData = 0:movieData.timeInterval_:(movieData.timeInterval_ * (nFrames-1));
    tUnits = 'Seconds';
else
    tData = 1:nFrames;
    tUnits = 'Frames';
end

%Figure handles for each movie
intFig = zeros(1,nChan);
bleachFig = zeros(1,nChan);
if hasMasks
    mIntFig = zeros(1,nChan);
    mBleachFig = zeros(1,nChan);
    mRatFig = zeros(1,nChan);
end

%% ------------------ Assessment ----------------------- %%

disp('Starting image assessment....')

for iChan = 1:nChan
    
    if ~p.BatchMode        
        waitbar((iChan-1)*nFrames / nImTot,wtBar,['Please wait, assessing channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end       
    
    % ------- Image Stats ------ %
    %Collect intensity statistics for each image
    
    disp(['Assessing images in ' imDirs{iChan} ])
    if hasMasks
        disp(['Using masks from ' maskDirs{iChan}]);
    end
    if hasBackMasks
        disp(['Using background masks from ' backMaskDirs{iChan}]);
    end
    
    for iFrame = 1:nFrames
        
        
        currImage = double(imread([imDirs{iChan} filesep imNames{iChan}{iFrame}]));
                
        imStats(iChan).AverageIntensity(iFrame) = mean(currImage(:));
        imStats(iChan).MinimumIntensity(iFrame) = min(currImage(:));
        imStats(iChan).MaximumIntensity(iFrame) = max(currImage(:));
        imStats(iChan).MedianIntensity(iFrame) = median(currImage(:));
        imStats(iChan).STDIntensity(iFrame) = std(currImage(:));
                
        
        if hasMasks
            currMask = imread([maskDirs{iChan} filesep maskNames{iChan}{iFrame}]);
            
            imStats(iChan).AverageMaskedIntensity(iFrame) = mean(currImage(currMask(:)));
            imStats(iChan).MinimumMaskedIntensity(iFrame) = min(currImage(currMask(:)));
            imStats(iChan).MaximumMaskedIntensity(iFrame) = max(currImage(currMask(:)));
            imStats(iChan).MedianMaskedIntensity(iFrame) = median(currImage(currMask(:)));
            imStats(iChan).STDMaskedIntensity(iFrame) = std(currImage(currMask(:)));
            imStats(iChan).AverageUnmaskedIntensity(iFrame) = mean(currImage(~currMask(:)));
            imStats(iChan).MinimumUnmaskedIntensity(iFrame) = min(currImage(~currMask(:)));
            imStats(iChan).MaximumUnmaskedIntensity(iFrame) = max(currImage(~currMask(:)));
            imStats(iChan).MedianUnmaskedIntensity(iFrame) = median(currImage(~currMask(:)));
            imStats(iChan).STDUnmaskedIntensity(iFrame) = std(currImage(~currMask(:)));            
            imStats(iChan).RatioMaskedToUnmaskedIntensity(iFrame) = imStats(iChan).AverageMaskedIntensity(iFrame) ...
                                                          / imStats(iChan).AverageUnmaskedIntensity(iFrame);
        
        end
        
        if hasBackMasks
            currBackMask = imread([backMaskDirs{iChan} filesep backMaskNames{iChan}{iFrame}]);
            
            imStats(iChan).AverageBackgroundIntensity(iFrame) = mean(currImage(currBackMask(:)));
            imStats(iChan).MinimumBackgroundIntensity(iFrame) = min(currImage(currBackMask(:)));
            imStats(iChan).MaximumBackgroundIntensity(iFrame) = max(currImage(currBackMask(:)));
            imStats(iChan).MedianBackgroundIntensity(iFrame) = median(currImage(currBackMask(:)));
            imStats(iChan).STDBackgroundIntensity(iFrame) = std(currImage(currBackMask(:)));           
            imStats(iChan).RatioMaskedToBackgroundIntensity(iFrame) = imStats(iChan).AverageMaskedIntensity(iFrame) ...
                                                          / imStats(iChan).AverageBackgroundIntensity(iFrame);
        end
        
        if ~p.BatchMode && mod(iFrame,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iFrame + (iChan-1)*nFrames) / nImTot,wtBar)
        end
        
        
    end
    
    % ------- Bleaching Calculations ------ %
    %Fit bleaching curves to data to get bleaching time constants
    
    %Fit function to intensity timeseries    
    [imStats(iChan).BleachingStats.FitCoef,resFit,~,covFit,mseFit] = nlinfit(tData(:),imStats(iChan).AverageIntensity(:),...
                                            fitFun,bInit,fitOptions);
    %Get confidence intervals of fit and fit values
    [imStats(iChan).BleachingStats.FitValues,imStats(iChan).BleachingStats.DeltaFit] = ...
                        nlpredci(fitFun,tData(:),imStats(iChan).BleachingStats.FitCoef,...
                        resFit,'covar',covFit,'mse',mseFit);
                    
    if hasMasks
        %Fit function to intensity timeseries
        [imStats(iChan).MaskedBleachingStats.FitCoef,resFit,~,covFit,mseFit] = ...
                    nlinfit(tData(:),imStats(iChan).AverageMaskedIntensity(:),...
                                                fitFun,bInit,fitOptions);
        %Get confidence intervals of fit and fit values
        [imStats(iChan).MaskedBleachingStats.FitValues,imStats(iChan).MaskedBleachingStats.DeltaFit] = ...
                            nlpredci(fitFun,tData(:),imStats(iChan).MaskedBleachingStats.FitCoef,...
                                resFit,'covar',covFit,'mse',mseFit);
    end
   
    
    
    % ------- Figure Making ------- %
    %Make, display and save figures for each channel
    
    
    %Whole-image intensity stats
	intFig(iChan) = figure(figOpts{:});
    hold on;
    errorbar(tData,imStats(iChan).AverageIntensity,imStats(iChan).STDIntensity);
    hold on
    plot(tData,imStats(iChan).MaximumIntensity,'--r');
    plot(tData,imStats(iChan).MinimumIntensity,'--g');
    plot(tData,imStats(iChan).MedianIntensity,'--m');
    legend('Average +- STD','Maximum','Minimum','Median')
    title(['Whole-Image Intensity Statistics, channel ' num2str(p.ChannelIndex(iChan))])
    xlabel(['Time, ' tUnits])
    ylabel('Intensity')
    xlim([min(tData) max(tData)])
    figName = [p.OutputDirectory filesep ...
        'whole image stats channel ' num2str(p.ChannelIndex(iChan))];
    hgsave(intFig(iChan),figName);
    print(intFig(iChan),figName,pOpt{:});
    
    %Whole-image bleaching stats
    bleachFig(iChan) = figure(figOpts{:});
    hold on;
    plot(tData,imStats(iChan).AverageIntensity,'k');
    plot(tData,imStats(iChan).BleachingStats.FitValues,'r')
    plot(tData,imStats(iChan).BleachingStats.FitValues + imStats(iChan).BleachingStats.DeltaFit,'--r')
    plot(tData,imStats(iChan).BleachingStats.FitValues - imStats(iChan).BleachingStats.DeltaFit,'--r')    
    legend('Average','Fit','Fit 95% C.I.')
    xlabel(['Time, ' tUnits ]);
    ylabel('Intensity')
    titleText{1} = ['Whole image bleaching statistics, channel ' num2str(p.ChannelIndex(iChan))];
    titleText{2} = ['y = ' num2str(imStats(iChan).BleachingStats.FitCoef(1))...
                ' * e^(' num2str(imStats(iChan).BleachingStats.FitCoef(2)) '*t)' ...
                ' + ' num2str(imStats(iChan).BleachingStats.FitCoef(3)) ' * ',...
                ' e^(' num2str(imStats(iChan).BleachingStats.FitCoef(4)) '*t)'];
    title(titleText,'Interpreter','none');
    figName = [p.OutputDirectory filesep ...
        'whole image bleaching stats channel ' num2str(p.ChannelIndex(iChan))];
    hgsave(bleachFig(iChan),figName);
    print(bleachFig(iChan),figName,pOpt{:});
    
    
    if hasMasks
        %Masked intensity stats
        mIntFig(iChan) = figure(figOpts{:});
        hold on;
        errorbar(tData,imStats(iChan).AverageMaskedIntensity,imStats(iChan).STDMaskedIntensity);
        hold on
        plot(tData,imStats(iChan).MaximumMaskedIntensity,'--r');
        plot(tData,imStats(iChan).MinimumMaskedIntensity,'--g');
        plot(tData,imStats(iChan).MedianMaskedIntensity,'--m');
        legend('Average +- STD','Maximum','Minimum','Median')
        title(['Masked Intensity Statistics, channel ' num2str(p.ChannelIndex(iChan))])
        xlabel(['Time, ' tUnits])
        ylabel('Masked Intensity')
        xlim([min(tData) max(tData)])
        figName = [p.OutputDirectory filesep ...
            'masked image stats channel ' num2str(p.ChannelIndex(iChan))];
        hgsave(mIntFig(iChan),figName);
        print(mIntFig(iChan),figName,pOpt{:});

        %Masked bleaching stats
        mBleachFig(iChan) = figure(figOpts{:});
        hold on;
        plot(tData,imStats(iChan).AverageMaskedIntensity,'k');
        plot(tData,imStats(iChan).MaskedBleachingStats.FitValues,'r')
        plot(tData,imStats(iChan).MaskedBleachingStats.FitValues + imStats(iChan).MaskedBleachingStats.DeltaFit,'--r')
        plot(tData,imStats(iChan).MaskedBleachingStats.FitValues - imStats(iChan).MaskedBleachingStats.DeltaFit,'--r')    
        legend('Average','Fit','Fit 95% C.I.')
        xlabel(['Time, ' tUnits ]);
        ylabel('Masked Intensity')
        titleText = cell(1,2);
        titleText{1} = ['Masked bleaching statistics, channel ' num2str(p.ChannelIndex(iChan))];
        titleText{2} = ['y = ' num2str(imStats(iChan).MaskedBleachingStats.FitCoef(1))...
                    ' * e^(' num2str(imStats(iChan).MaskedBleachingStats.FitCoef(2)) '*t)' ...
                    ' + ' num2str(imStats(iChan).MaskedBleachingStats.FitCoef(3)) ' * ',...
                    ' e^(' num2str(imStats(iChan).MaskedBleachingStats.FitCoef(4)) '*t)'];
        title(titleText,'Interpreter','none');
        figName = [p.OutputDirectory filesep ...
            'masked image bleaching stats channel ' num2str(p.ChannelIndex(iChan))];
        hgsave(mBleachFig(iChan),figName);
        print(mBleachFig(iChan),figName,pOpt{:});
        
        %Masked/unmasked ratio figure
        mRatFig(iChan) = figure(figOpts{:});
        hold on;
        plot(tData,imStats(iChan).RatioMaskedToUnmaskedIntensity)
        legNames = {'Masked:Unmasked'};
        if hasBackMasks
            plot(tData,imStats(iChan).RatioMaskedToBackgroundIntensity,'r')
            legNames = [legNames {'Masked:Background'}]; %#ok<AGROW>
        end
        xlabel(['Time, ' tUnits])
        ylabel('Intensity Ratio')
        legend(legNames)
        title(['Ratios of average segmented intensities, channel ' num2str(p.ChannelIndex(iChan))])
        figName = [p.OutputDirectory filesep ...
            'masked image ratios channel ' num2str(p.ChannelIndex(iChan))];
        hgsave(mRatFig(iChan),figName);
        print(mRatFig(iChan),figName,pOpt{:});                           
        
    end
    
               
    
end


%% ------------------- Finalize ------------------ %%


%Close all figures if in batch mode. Otherwise leave them open for user to
%look at.
if p.BatchMode
    for j = 1:nChan
        close(intFig(j))
        close(bleachFig(j))
        if hasMasks
            close(mIntFig(iChan))
            close(mBleachFig(iChan))
            close(mRatFig(iChan))
        end
    end
end

%Save the calculated statistics for each channel
for j = 1:nChan
    %Get the current channel because the save function is weird
    imageStats = imStats(iChan); %#ok<NASGU>
    save([p.OutputDirectory filesep 'image statistics channel ' num2str(p.ChannelIndex(j))],'imageStats');    
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end
disp('Finished image assessment....')






