function [data] = determinePitDensities(data)
% determinePitDensities calculates density of pits in whole movie
% 
% INPUT:    data  = struct containing all data for one condition (e.g. 
%                       clathrin control)
%           
%
% OUTPUT:   data.pitDensity [av std]
%
% Dinah Loerke, last changed 03/19/2008

od = cd;

for i = 1:length(data)
    
    currDir = data(i).source;
    cd(currDir);
    
    if exist('TrackInfoMatrices');
        
        cd('TrackInfoMatrices');
        ti = load('trackInfo.mat');
        if isfield(ti,'trackInfo')
            currTrackInfo = ti.trackInfo;
        elseif isfield(ti,'trackInfoMat')
            currTrackInfo = ti.trackInfoMat;
        end
        
        [stx,sty] = size(currTrackInfo);
        numf = sty/8;
        
        % extract all x and y positions of object
        cmx = full(currTrackInfo(:,1:8:sty));
        cmy = full(currTrackInfo(:,2:8:sty));
        
        % number of objects detected in each frame
        for k=1:numf
            numPitsFrame(k) = length(find(cmx(:,k)>0));
        end
        
        % all positions together
        allx = nonzeros(cmx(:));
        ally = nonzeros(cmy(:));
        nx = length(allx(:));
        ny = length(ally(:));
        numPits = nx/numf;
        
        % earlier version used convex hull; this can be tricky for more
        % complicated cell shaped, so allow more complex edge
        % calculate 'basis mask', which is the overlay of all detected objects
        % first, set central pixels of detected objects to 1
        mask_detection = zeros( round(max(allx)), round(max(ally))  );
        for t=1:length(allx)
            detxpos = max(1,round(allx(t)));
            detypos = max(1,round(ally(t)));
            mask_detection(detxpos,detypos) = 1;
        end
        % second, dilate the single pixels to little discs
        se5 = strel('disk',5);
        se20 = strel('disk',20);
        mask_detectdilate = imdilate(logical(mask_detection),se5);
        % third, fill the resulting area and filter with Gaussian profile to
        % create a continuous area
        mask_filled = imfill(double(mask_detectdilate),'holes');
        mask_close = imclose(mask_filled,se20);
        
        totalarea = sum(mask_close(:));
        
        denPitsFrame = 1000*numPitsFrame/totalarea;        
              
        data(i).pitDensity = denPitsFrame;
        localdensity(i) = nanmean(denPitsFrame);
    else
        data(i).pitDensity = nan;
    end
    
    cd(od);
    
end % of for m

avDensity = round(100*nanmean(localdensity))/100;
stdDensity = round(100*nanstd(localdensity))/100;
disp(['pit density = ',num2str(avDensity),' +- ',num2str(stdDensity)]);          
        
end % of function
    
   
    
    