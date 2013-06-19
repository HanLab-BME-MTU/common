function plotPropertySpatialMap2D(tracksFinal,diffAnalysisRes,diffModeAnRes,...
    numNeighbors,properties2plot,positions2plot,lengthMinMax,fixedSegLength,...
    figureName,image)
%PLOTPROPERTYSPATIALMAP creates spatial maps of trajectory properties
%
%SYNOPSIS plotPropertySpatialMap2D(tracksFinal,diffAnalysisRes,diffModeAnRes,...
%    numNeighbors,properties2plot,positions2plot,lengthMinMax,fixedSegLength,...
%    figureName,image)
%
%INPUT  tracksFinal    : Output of trackCloseGapsKalman.
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%       diffModeAnRes  : Output of trackDiffModeAnalysis.
%       numNeighbors   : Output of numNeighborsTrack.
%       properties2plot: Row vector of properties to plot:
%                        1 - Diffusion classification.
%                        2 - Diffusion coefficient from MSS analysis.
%                        3 - Confinement radius.
%                        4 - Track lifetime.
%                        5 - Average frame-to-frame displacement.
%                        6 - Diffusion mode.
%                        7 - Diffusion coefficient from mode analysis.
%                        8 - Number of neighbors.
%                        Optional. Default: all.
%	    positions2plot : Row vector of trajectory position to plot:
%                        1 - Center position.
%                        2 - Start position.
%                        3 - End position.
%                        Optional. Default: 1.
%       lengthMinMax   : Minimum and maximum length of trajectories to
%                        include in plots.
%                        Optional. Default: [5 99].
%       fixedSegLength : 1 to divide a property range into segments of
%                        fixed length, 0 to divide a property range into
%                        segments of variable length such that they contain
%                        an equal number of elements.
%                        Optional. Default: 1;
%       figureName     : Figure name. 
%                        Optional. Default: [].
%       image          : Image(s) to overlay spatial map on.
%                        Can be 1 or 2 images.
%                        If 1 image, the image will be displayed as is.
%                        If 2 images (usually masks), then only their
%                        outline will be plotted (green for image 1, red
%                        for image 2)
%                        Optional. Default: None.
%
%REMARKS If there are n properties2plot and m positions2plot, the code will
%output n*m spatial maps.
%
%Khuloud Jaqaman, August 2009

%% Input

if nargin < 4
    disp('--plotPropertySpatialMap2D: Incorrect number of input arguments!');
    return
end

if nargin < 5 || isempty(properties2plot)
    properties2plot = 1:8;
end

if nargin < 6 || isempty(positions2plot)
    positions2plot = 1;
end

if nargin < 7 || isempty(lengthMinMax)
    lengthMinMax = [5 99];
end

if nargin < 8 || isempty(fixedSegLength)
    fixedSegLength = 1;
end

if nargin < 9 || isempty(figureName)
    figureName = [];
end

if nargin < 10 || isempty(image)
    image = [];
end

%% Preparation for plotting

%define the default of 10 segments for plotting some of the properties
numSegments = 10;

%assign segment colors
segmentColor = [0 0 0; 0 0 1; 0.2 0.7 0.7; 0 1 1; 0 1 0; ...
    0.6824 0.4667 0; 1 0.7 0; 1 0 0; 1 0 1; 0.7 0 1];

%construct parts of figure titles
plottedProperty = {'MSS classification','MSS diffusion coefficient',...
    'Confinement radius','Lifetime','Frame-to-frame displacement',...
    'Diffusion mode','Mode analysis diffusion coefficient','Number of neighbors'};
plottedPosition = {'center position','start position','end position'};

%% Trajectory pre-processing

%keep only trajectories of acceptable length
criteria.lifeTime.min = lengthMinMax(1);
criteria.lifeTime.max = lengthMinMax(2);
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx);
diffAnalysisRes = diffAnalysisRes(indx);
diffModeAnRes = diffModeAnRes(indx);
numNeighbors = numNeighbors(indx);

%save tracksFinal into a new variable name
inputStructure = tracksFinal;

%if tracksFinal is a structure ...
if isstruct(tracksFinal)
    
    %get number of compound tracks
    numCompTrack = length(tracksFinal);
    
    %get the number of segments per compound track
    numTrackSeg = getNumSegmentsPerTrack(tracksFinal);
    
    %determine the location of the 1st segment of each compound track in
    %the "big matrix" of all segments
    segLoc = [1; cumsum(numTrackSeg(1:end-1))+1];
    
    %get the start, end and life time information of each segment
    trajSEL = getTrackSEL(tracksFinal,1);
    
    %from this get total number of segments (called "traj")
    numTraj = size(trajSEL,1);
    
    %get the start, end and center position of trajectories (segments)
    %also calculate the average frame-to-frame displacement
    trajXCoord = NaN(numTraj,3);
    trajYCoord = NaN(numTraj,3);
    frame2frameDisp = NaN(numTraj,1);
    for iCompTrack = 1 : numCompTrack
        
        %get current compound track's information over its lifetime
        infoCompTrack = tracksFinal(iCompTrack).tracksCoordAmpCG;
        xCoordCompTrack = infoCompTrack(:,1:8:end);
        yCoordCompTrack = infoCompTrack(:,2:8:end);
        
        %get "local" start, end and life times
        segSEL = getTrackSEL(infoCompTrack);
        
        for iSegment = 1 : numTrackSeg(iCompTrack)
            
            %get current segment's positions over its lifetime
            xCoordCurrent = xCoordCompTrack(iSegment,segSEL(iSegment,1):segSEL(iSegment,2));
            yCoordCurrent = yCoordCompTrack(iSegment,segSEL(iSegment,1):segSEL(iSegment,2));
            
            %calculate start, end and center positions
            startPos = [xCoordCurrent(1) yCoordCurrent(1)];
            endPos = [xCoordCurrent(end) yCoordCurrent(end)];
            centerPos = [nanmean(xCoordCurrent) nanmean(yCoordCurrent)];
            
            %assemble the position information, ordered as instructed for the input
            %variable positions2plot
            trajXCoord(segLoc(iCompTrack)+iSegment-1,:) = [centerPos(1) startPos(1) endPos(1)];
            trajYCoord(segLoc(iCompTrack)+iSegment-1,:) = [centerPos(2) startPos(2) endPos(2)];
            
            %calculate the average frame-to-frame displacement (root mean
            %square)
            frame2frameDisp(segLoc(iCompTrack)+iSegment-1) = sqrt(nanmean(diff(xCoordCurrent).^2+diff(yCoordCurrent).^2));
            
        end
    end
    
else %if tracksFinal is a matrix ...
    
    %get number of trajectories
    numTraj = size(tracksFinal,1);
    
    %extract the x- and y-coordinates from the big matrix
    xCoord = tracksFinal(:,1:8:end);
    yCoord = tracksFinal(:,2:8:end);
    
    %get the start, end and life time information of trajectories
    trajSEL = getTrackSEL(tracksFinal);
    
    %get the start, end and center position of trajectories
    %also calculate the average frame-to-frame displacement
    trajXCoord = NaN(numTraj,3);
    trajYCoord = NaN(numTraj,3);
    frame2frameDisp = NaN(numTraj,1);
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
        
        %calculate the average frame-to-frame displacement (root mean
        %square)
        frame2frameDisp(iTraj) = sqrt(nanmean(diff(xCoordCurrent).^2+diff(yCoordCurrent).^2));
        
    end
    
end

%find x-coordinate limits
minXCoord = min(floor(min(trajXCoord(:))),0);
maxXCoord =  ceil(max(trajXCoord(:)));

%find y-coordinate limits
minYCoord = min(floor(min(trajYCoord(:))),0);
maxYCoord =  ceil(max(trajYCoord(:)));

%% Property extraction and pre-processing

if any(properties2plot==1)
    
    %get classifications from MSS diffusion analysis results
    trajClass = vertcat(diffAnalysisRes.classification);
    trajClass = trajClass(:,2);
    
    %calculate the fraction of trajectories in each classification
    fracTrajClass = hist(trajClass,1:3);
    fracTrajClass = fracTrajClass / sum(fracTrajClass);
    
end

if any(properties2plot==2)
    
    %get diffusion coefficients from diffusion analysis results
    % diffCoefNorm = catStruct(1,'diffAnalysisRes.fullDim.normDiffCoef');
    diffCoefGen = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');
    
    %divide the range of diffusion coefficients into segments, determine
    %which segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    % [diffCoefSegment,segmentEdgesDC,fracInSegmentsDC] = divideRangeIntoSegments(...
    %     diffCoefNorm,numSegments,fixedSegLength);
    [diffCoefSegment,segmentEdgesDC,fracInSegmentsDC] = divideRangeIntoSegments(...
        diffCoefGen,numSegments,fixedSegLength);
    
end

if any(properties2plot==3)
    
    %get confinement radii from diffusion analysis results
    confRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');
    
    %divide the range of confinement radii into segments, determine which
    %segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [confRadSegment,segmentEdgesCR,fracInSegmentsCR] = divideRangeIntoSegments(...
        confRad,numSegments,fixedSegLength);
    
end

if any(properties2plot==4)
    
    %get trajectory lifetimes
    trajLft = trajSEL(:,3);
    
    %divide the range of lifetimes into segments, determine which
    %segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [trajLftSegment,segmentEdgesLft,fracInSegmentsLft] = divideRangeIntoSegments(...
        trajLft,numSegments,fixedSegLength);
    
end

if any(properties2plot==5)
    
    %divide the range of frame-to-frame displacements into segments, determine
    %which segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [f2fDispSegment,segmentEdgesFFD,fracInSegmentsFFD] = divideRangeIntoSegments(...
        frame2frameDisp,numSegments,fixedSegLength);
    
end

if any(properties2plot==6)
    
    %get diffusion modes from diffusion mode analysis results
    trajDiffMode = vertcat(diffModeAnRes.diffMode);
    
    %calculate the fraction of trajectories in each mode
    fracDiffMode = hist(trajDiffMode,1:4);
    fracDiffMode = fracDiffMode / sum(fracDiffMode);
    
end

if any(properties2plot==7)
    
    %get diffusion coefficients from diffusion mode analysis results
    diffCoefDisp2 = vertcat(diffModeAnRes.diffCoef);
    
    %divide the range of diffusion coefficients into segments, determine
    %which segment each trajectory falls into, and calculate the fraction of
    %trajectories in each segment
    [diffCoef2Segment,segmentEdgesDC2,fracInSegmentsDC2] = divideRangeIntoSegments(...
        diffCoefDisp2,numSegments,fixedSegLength);
    
end

if any(properties2plot==8)
    
    %get number of neighbors from numNeighbors variable
    smDensity = vertcat(numNeighbors.value);
    
    %divide the range of number of neighbors into segements, determine which segment
    %each trajectory falls into, and calculate the fraction of trajectories
    %in each segement
    [smDensitySegment,segmentEdgesDensity,fracInSegmentsDensity] = divideRangeIntoSegments(...
        smDensity,numSegments,fixedSegLength);
    
end

% %% Some numbers for output
% 
% %get indices of trajectories in the various classifications
% indxConf = find(trajClass==1);
% indxBrown = find(trajClass==2);
% indxDir = find(trajClass==3);
% 
% f2fDispRMS = [nanmean(frame2frameDisp(indxConf)) nanmean(frame2frameDisp(indxBrown)) ...
%     nanmean(frame2frameDisp(indxDir)) nanmean(frame2frameDisp([indxConf;indxBrown;indxDir])) ...
%     nanmean(frame2frameDisp)];

%% Plotting

%go over all properties and positions
for iProperty = properties2plot
    for iPos = positions2plot
        
        %make new figure
        if isempty(figureName)
            h = figure;
        else
            h = figure('Name',figureName);
        end

        % SUBPLOT 1: Spatial map of values %
        
        subplot(1,3,[1 2])
        hold on
        
        %plot the image to overlay the spatial map on, if given
        if size(image,3)==1
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
                                
            case 2 %diffusion coefficient MSS
                
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
                    end
                end
                
            case 4 %lifetime
                
                %go over the different lifetime segments
                for iSegment = 1 : numSegments
                    indx = find(trajLftSegment==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                            '.','Color',segmentColor(iSegment,:));
                    end
                end
                
            case 5 %frame-to-frame displacement
                
                %trajectories without a frame-to-frame displacement
                indx = find(isnan(f2fDispSegment));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end
                
                %go over the different frame-to-frame displacement segments
                for iSegment = 1 : numSegments
                    indx = find(f2fDispSegment==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                            '.','Color',segmentColor(iSegment,:));
                    end
                end
                
            case 6 %diffusion mode
                
                %unclassified
                indx = find(isnan(trajDiffMode));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end
                
                %mode 3
                indx = find(trajDiffMode == 3);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','c');
                end
                
                %mode 2
                indx = find(trajDiffMode == 2);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','b');
                end
                
                %mode 4
                indx = find(trajDiffMode == 4);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','m');
                end
                
                %mode 1
                indx = find(trajDiffMode == 1);
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color','k');
                end
                
            case 7 %diffusion coefficient from mode analysis
                
                %trajectories without a diffusion coefficient i.e.
                %unclassified trajectories
                indx = find(isnan(diffCoef2Segment));
                if ~isempty(indx)
                    plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                        '.','Color',[0.7 0.7 0.7]);
                    legendText{end+1} = 'unclassified'; %#ok<AGROW>
                end
                
                %go over the different diff. coef. segments
                for iSegment = 1 : numSegments
                    indx = find(diffCoef2Segment==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                            '.','Color',segmentColor(iSegment,:));
                    end
                end
                
            case 8 %number of neighbors
                
                %go over the different density segments
                for iSegment = 1 : numSegments
                    indx = find(smDensitySegment==iSegment);
                    if ~isempty(indx)
                        plot(trajXCoord(indx,iPos),trajYCoord(indx,iPos),...
                            '.','Color',segmentColor(iSegment,:));
                    end
                end
                
        end
        
        if size(image,3)==2
            maskBounds = bwboundaries(image(:,:,1));
            cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',3)),maskBounds);
            maskBounds = bwboundaries(image(:,:,2));
            cellfun(@(x)(plot(x(:,2),x(:,1),'r--','LineWidth',3)),maskBounds);
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
        ylabel('Fraction');
        
        %initialize figure legend text
        legendText = [];
        
        switch iProperty
            
            case 1 %classification
                
                %plot the bars with different colors
                bar([1 2],[fracTrajClass(1) 0],'BarWidth',1,'FaceColor',...
                    'b','EdgeColor','none');
                bar([2 3],[fracTrajClass(2) 0],'BarWidth',1,'FaceColor',...
                    'c','EdgeColor','none');
                bar([3 4],[fracTrajClass(3) 0],'BarWidth',1,'FaceColor',...
                    'm','EdgeColor','none');
                legendText = {'confined','free','directed'};
                
            case 2 %diffusion coefficient MSS
                
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
                
            case 3 %confinement radius
                
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
                
            case 4 %lifetime
                
                %get the center and width of each segment
                segmentCenter = mean(segmentEdgesLft,2);
                segmentWidth = segmentEdgesLft(:,2) - segmentEdgesLft(:,1);
                
                %plot the bars with different colors
                for iSegment = 1 : numSegments
                    bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
                        segmentWidth(iSegment)],[fracInSegmentsLft(iSegment) 0],...
                        'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
                        'EdgeColor','none');
                    legendText{end+1} = [num2str(segmentEdgesLft(iSegment,1)) ...
                        ' - ' num2str(segmentEdgesLft(iSegment,2))]; %#ok<AGROW>
                end
                
            case 5 %frame-to-frame displacement
                
                %get the center and width of each segment
                segmentCenter = mean(segmentEdgesFFD,2);
                segmentWidth = segmentEdgesFFD(:,2) - segmentEdgesFFD(:,1);
                
                %plot the bars with different colors
                for iSegment = 1 : numSegments
                    bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
                        segmentWidth(iSegment)],[fracInSegmentsFFD(iSegment) 0],...
                        'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
                        'EdgeColor','none');
                    legendText{end+1} = [num2str(segmentEdgesFFD(iSegment,1)) ...
                        ' - ' num2str(segmentEdgesFFD(iSegment,2))]; %#ok<AGROW>
                end
                
            case 6 %diffusion mode
                
                %plot the bars with different colors
                bar([1 2],[fracDiffMode(1) 0],'BarWidth',1,'FaceColor',...
                    'k','EdgeColor','none');
                bar([2 3],[fracDiffMode(2) 0],'BarWidth',1,'FaceColor',...
                    'b','EdgeColor','none');
                bar([3 4],[fracDiffMode(3) 0],'BarWidth',1,'FaceColor',...
                    'c','EdgeColor','none');
                bar([4 5],[fracDiffMode(4) 0],'BarWidth',1,'FaceColor',...
                    'm','EdgeColor','none');
                legendText = {'mode 1','mode 2','mode 3','mode 4'};
                
            case 7 %diffusion coefficient from mode analysis
                
                %get the center and width of each segment
                segmentCenter = mean(segmentEdgesDC2,2);
                segmentWidth = segmentEdgesDC2(:,2) - segmentEdgesDC2(:,1);
                
                %plot the bars with different colors
                for iSegment = 1 : numSegments
                    bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
                        segmentWidth(iSegment)],[fracInSegmentsDC2(iSegment) 0],...
                        'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
                        'EdgeColor','none');
                    legendText{end+1} = [num2str(segmentEdgesDC2(iSegment,1)) ...
                        ' - ' num2str(segmentEdgesDC2(iSegment,2))]; %#ok<AGROW>
                end
                
            case 8 %number of neighbors
                
                %get the center and width of each segment
                segmentCenter = mean(segmentEdgesDensity,2);
                segmentWidth = segmentEdgesDensity(:,2) - segmentEdgesDensity(:,1);
                
                %plot the bars with different colors
                for iSegment = 1 : numSegments
                    bar([segmentCenter(iSegment) segmentCenter(iSegment)+...
                        segmentWidth(iSegment)],[fracInSegmentsDensity(iSegment) 0],...
                        'BarWidth',1,'FaceColor',segmentColor(iSegment,:),...
                        'EdgeColor','none');
                    legendText{end+1} = [num2str(segmentEdgesDensity(iSegment,1)) ...
                        ' - ' num2str(segmentEdgesDensity(iSegment,2))]; %#ok<AGROW>
                end
                
        end
        
        %add subplot title
        title(['Distribution of ' plottedProperty{iProperty} ' values'])
        
        %write color legend - this won't include unclassified/unconfined
        legend(legendText)
        
        %hold off
        hold off
        
        %save figure
        saveas(h,[figureName '_' plottedProperty{iProperty} '_vs_' plottedPosition{iPos}],'fig');
        
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

