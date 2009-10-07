function plotPropertySpatialMap2D(tracksFinal,figureName,minLength,...
    properties2plot,positions2plot,image,fixedSegLength,diffAnalysisRes)
%PLOTPROPERTYSPATIALMAP creates spatial maps of trajectory properties
%
%SYNOPSIS plotPropertySpatialMap2D(tracksFinal,figureName,minLength,...
%    properties2plot,positions2plot,image,fixedSegLength,diffAnalysisRes)
%
%INPUT  tracksFinal    : Output of trackCloseGapsKalman.
%       figureName     : Figure name.
%       minLength      : Minimum length of a trajectory to analyze and
%                        include in plots.
%                        Optional. Default: 20.
%       properties2plot: Row vector of properties to plot:
%                        1 - Diffusion classification.
%                        2 - Diffusion coefficient.
%                        3 - Confinement radius.
%                        Optional. Default: Everything.
%	    positions2plot : Row vector of trajectory position to plot:
%                        1 - Center position.
%                        2 - Start position.
%                        3 - End position.
%                        Optional. Default: 1.
%       image          : Image to overlay spatial map on.
%                        Optional. Default: [].
%       fixedSegLength : 1 to divide a property range into segments of
%                        fixed length, 0 to divide a property range into
%                        segments of variable length such that they contain
%                        an equal number of elements.
%                        Optional. Default: 1;
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%                        Optional. If not input but needed, it will be
%                        calculated within the code.
%
%REMARKS If there are n properties2plot and m positions2plot, the code will
%output n*m spatial maps.
%
%Khuloud Jaqaman, August 2009

%% Input

if nargin < 1
    disp('--plotPropertySpatialMap2D: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(figureName)
    figureName = [];
end

if nargin < 3 || isempty(minLength)
    minLength = 20;
end

if nargin < 4 || isempty(properties2plot)
    properties2plot = 1:3;
end

if nargin < 5 || isempty(positions2plot)
    positions2plot = 1;
end

if nargin < 6 || isempty(image)
    image = [];
end

if nargin < 7 || isempty(fixedSegLength)
    fixedSegLength = 1;
end

if nargin < 8 || isempty(diffAnalysisRes)
    diffAnalysisRes = [];
end

%% Preparation for plotting

%define the default of 10 segments for plotting some of the properties
numSegments = 10;

%assign segment colors
segmentColor = [0 0 0; 0 0 1; 0.2 0.7 0.7; 0 1 1; 0 1 0; ...
    0.6824 0.4667 0; 1 0.7 0; 1 0 0; 1 0 1; 0.7 0 1];

%construct parts of figure titles
plottedProperty = {'Classification','Diffusion coefficient',...
    'Confinement radius'};
plottedPosition = {'center position','start position','end position'};

%% Trajectory pre-processing

%keep only trajectories longer than minLength
criteria.lifeTime.min = minLength;
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx);

%convert tracksFinal into matrix if it's a structure
inputStructure = tracksFinal;
if isstruct(tracksFinal)
    clear tracksFinal
    tracksFinal = convStruct2MatIgnoreMS(inputStructure);
end

%get number of trajectories
numTraj = size(tracksFinal,1);

%extract the x- and y-coordinates from the big matrix
xCoord = tracksFinal(:,1:8:end);
yCoord = tracksFinal(:,2:8:end);

%find x-coordinate limits
minXCoord = min(floor(min(xCoord(:))),0);
maxXCoord =  ceil(max(xCoord(:)));

%find y-coordinate limits
minYCoord = min(floor(min(yCoord(:))),0);
maxYCoord =  ceil(max(yCoord(:)));

%get the start, end and life time information of trajectories
trajSEL = getTrackSEL(tracksFinal);

%get the start, end and center position of trajectories
trajXCoord = NaN(numTraj,3);
trajYCoord = NaN(numTraj,3);
for iTraj = 1 : numTraj

    %get current track's positions over its lifetime
    xCoordCurrent = xCoord(iTraj,trajSEL(iTraj,1):trajSEL(iTraj,2));
    yCoordCurrent = yCoord(iTraj,trajSEL(iTraj,1):trajSEL(iTraj,2));

    %calculate start, end and center positions
    startPos = [xCoordCurrent(1) yCoordCurrent(1)];
    endPos = [xCoordCurrent(end) yCoordCurrent(end)];
    centerPos = [nanmean(xCoordCurrent) nanmean(yCoordCurrent)];

    %assemble the position information, ordered as instructed for the input
    %variable positions2plot
    trajXCoord(iTraj,:) = [centerPos(1) startPos(1) endPos(1)];
    trajYCoord(iTraj,:) = [centerPos(2) startPos(2) endPos(2)];

end

%% Property extraction and pre-processing

%perform diffusion analysis if not supplied
if any(properties2plot<=3) && isempty(diffAnalysisRes)
    diffAnalysisRes = trackDiffusionAnalysis1(inputStructure,1,2,0,0.05);
end

%get classifications from diffusion analysis results
trajClass = vertcat(diffAnalysisRes.classification);
trajClass = trajClass(:,2);

%calculate the fraction of trajectories in each classification
fracTrajClass = hist(trajClass,1:3);
fracTrajClass = fracTrajClass / sum(fracTrajClass);

%get diffusion coefficients from diffusion analysis results
diffCoefNorm = catStruct(1,'diffAnalysisRes.fullDim.normDiffCoef');

%divide the range of diffusion coefficients into segments, determine
%which segment each trajectory falls into, and calculate the fraction of
%trajectories in each segment
[diffCoefSegment,segmentEdgesDC,fracInSegmentsDC] = divideRangeIntoSegments(...
    diffCoefNorm,numSegments,fixedSegLength);

%get confinement radii from diffusion analysis results
confRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');

%divide the range of confinement radii into segments, determine which
%segment each trajectory falls into, and calculate the fraction of
%trajectories in each segment
[confRadSegment,segmentEdgesCR,fracInSegmentsCR] = divideRangeIntoSegments(...
    confRad,numSegments,fixedSegLength);

%% Plotting

%go over all properties and positions
for iProperty = properties2plot
    for iPos = positions2plot
        
        %make new figure
        if isempty(figureName)
            figure
        else
            figure('Name',figureName)
        end

        % SUBPLOT 1: Spatial map of values %
        
        subplot(1,3,[1 2])
        hold on
        
        %plot the image to overlay the spatial map on, if given
        if ~isempty(image)
            imshow(image,[]);
        else
            imshow(ones(maxYCoord,maxXCoord),[]);
        end
        
        %set figure axes limits
        axis([minXCoord maxXCoord minYCoord maxYCoord]);
        
        %show coordinates on axes
        axH = gca;
        set(axH,'visible','on');
        
        %label axes
        xlabel('x-coordinate (pixels)');
        ylabel('y-coordinate (pixels)');
        
        %initialize figure legend text
        legendText = [];
        
        %make plots
        switch iProperty
            
            case 1 %classification
                
                %unclassified
                indx = find(isnan(trajClass));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end

                %free diffusion
                indx = find(trajClass == 2);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','c');
                end
                
                %confined
                indx = find(trajClass == 1);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','b');
                end
                
                %directed
                indx = find(trajClass == 3);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','m');
                end
                                
            case 2 %diffusion coefficient
                
                %trajectories without a diffusion coefficient i.e.
                %unclassified trajectories
                indx = find(isnan(diffCoefSegment));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end
                
                %go over the different diff. coef. segments
                for iSegment = 1 : numSegments
                    indx = find(diffCoefSegment==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                            '.','Color',segmentColor(iSegment,:));
                        %                         legendText{end+1} = [num2str(segmentEdgesDC(iSegment,1)) ...
                        %                             ' - ' num2str(segmentEdgesDC(iSegment,2)) '; ' num2str(fracInSegmentsDC(iSegment))];
                    end
                end
                
            case 3 %confinement radius
                
                %unclassified trajectories
                indx = find(isnan(confRadSegment) & isnan(trajClass));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end
                
                %trajectories that are not confined
                indx = find(isnan(confRadSegment) & ~isnan(trajClass));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '+','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'not confined'; %#ok<AGROW>
                end
                
                %go over the different conf. rad. segment
                for iSegment = 1 : numSegments
                    indx = find(confRadSegment==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                            '.','Color',segmentColor(iSegment,:));
                        %                         legendText{end+1} = [num2str(segmentEdgesCR(iSegment,1)) ...
                        %                             ' - ' num2str(segmentEdgesCR(iSegment,2))];
                    end
                end
                
        end
        
        %add subplot title
        title([plottedProperty{iProperty} ' vs. ' plottedPosition{iPos}])
        
        %write color legend - this will be only for unclassified (and
        %unconfined in the case of confinement radius)
        legend(legendText)
        
        %hold off
        hold off
        
        % SUBPLOT 2: Distribution of values %
        
        subplot(1,3,3)
        hold on
        
        %label axes
        xlabel(plottedProperty{iProperty});
        ylabel('Fraction of trajectories');
        
        %initialize figure legend text
        legendText = [];
        
        switch iProperty
            
            case 1
                
                %plot the bars with different colors
                bar([1 2],[fracTrajClass(1) 0],'BarWidth',1,'FaceColor',...
                    'b','EdgeColor','none');
                bar([2 3],[fracTrajClass(2) 0],'BarWidth',1,'FaceColor',...
                    'c','EdgeColor','none');
                bar([3 4],[fracTrajClass(3) 0],'BarWidth',1,'FaceColor',...
                    'm','EdgeColor','none');
                legendText = {'confined','free','directed'};
                
            case 2
                
                %get the center and width of each segment
                segmentCenter = mean(segmentEdgesDC,2);
                segmentWidth = segmentEdgesDC(:,2) - segmentEdgesDC(:,1);
                
                %plot the bars with different colors
                for iSegment = 1 : numSegments
                    bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
                        segmentWidth(iSegment)],[fracInSegmentsDC(iSegment) 0],...
                        'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
                        'EdgeColor','none');
                    legendText{end+1} = [num2str(segmentEdgesDC(iSegment,1)) ...
                        ' - ' num2str(segmentEdgesDC(iSegment,2))]; %#ok<AGROW>
                end
                
            case 3
                
                %get the center and width of each segment
                segmentCenter = mean(segmentEdgesCR,2);
                segmentWidth = segmentEdgesCR(:,2) - segmentEdgesCR(:,1);
                
                %plot the bars with different colors
                for iSegment = 1 : numSegments
                    bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
                        segmentWidth(iSegment)],[fracInSegmentsCR(iSegment) 0],...
                        'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
                        'EdgeColor','none');
                    legendText{end+1} = [num2str(segmentEdgesCR(iSegment,1)) ...
                        ' - ' num2str(segmentEdgesCR(iSegment,2))]; %#ok<AGROW>
                end
                
        end
        
        %add subplot title
        title(['Distribution of ' plottedProperty{iProperty} ' values'])
        
        %write color legend - this won't include unclassified/unconfined
        legend(legendText)
        
        %hold off
        hold off
        
    end
end




%% subfunctions

function [valVectorSegment,segmentEdges,fracInSegment] = divideRangeIntoSegments(...
    valVector,numSegments,fixedSegLength)

%find minimum and maximum values
minVal = min(valVector);
maxVal = max(valVector);

if fixedSegLength %if fixed segment length ...
    
    %divide the range of values into equal segments
    valIncrement = (maxVal - minVal) / numSegments;
    
    %define the segment upper edges
    segmentUpperEdges = minVal + (1 : numSegments) * valIncrement;
    
else %if equal distribution of elements among segments ...
    
    %determine segment upper edges based on percentiles
    segmentUpperEdges = prctile(valVector,(1:numSegments)*100/numSegments);
    
end

%determine which segment each value falls into
valVectorSegment = NaN(size(valVector));
for iSegment = numSegments : -1 : 1
    valVectorSegment(valVector<=segmentUpperEdges(iSegment)) = iSegment;
end

%store the segment edges
segmentEdges = [[minVal; segmentUpperEdges(1:end-1)'] segmentUpperEdges'];

%calculate the fraction of values in each segment
fracInSegment = hist(valVectorSegment,1:numSegments);
fracInSegment = fracInSegment/sum(fracInSegment);


% function segmentColor = distributeSegmentColors(numSegments)
% 
% %define the color regimes - fixed to 5 for now
% numColors = 5;
% regimeEdgeColor = [0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1]; %bcgyrm
% 
% %generate vector of segment indices
% segmentIndx = 1 : numSegments;
% 
% %initialize the vector of segment colors
% segmentColor = zeros(numSegments,3);
% 
% %go over color regimes and assign segments to each color regime
% for iColor = 1 : numColors
% 
%     %divide number of remaining segments by number of remaining colors
%     numSegments4thisColor = length(segmentIndx) / (numColors - iColor + 1);
%     
%     %assign segments to this color regime
%     segments4thisColor = segmentIndx(1:numSegments4thisColor);
%     
%     %assign colors to the segments of this color regime
%     for iSegment = 1 : numSegments4thisColor
%         segmentColor(segments4thisColor(iSegment),:) = ...
%             regimeEdgeColor(iColor,:) * ...
%             (numSegments4thisColor-iSegment+1)/numSegments4thisColor ...
%             + regimeEdgeColor(iColor+1,:) * (iSegment-1)/numSegments4thisColor;
%     end
% 
%     %update vector of segment indices by removing the indices already
%     %assigned to this color
%     segmentIndx = segmentIndx(numSegments4thisColor+1:end);
%     
% end

