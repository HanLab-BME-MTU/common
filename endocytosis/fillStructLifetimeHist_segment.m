function [data] = fillStructLifetimeHist_segment(data);
% fill experiment structure with lifetime histogram based on the lftInfo
% file saved at the specified directory; separate lftHist vectors are
% calculated based on the image segmentation in the specified folders
% 
% SYNOPSIS [exp]=fillStructLifetimeInfo(exp);
%
% INPUT     data: experiment structure, which has to contain the fields
%                .source
%                .movieLength
%                .segmentDataFileName 
%                .segmentDataFilePath 
% OUTPUT    creates new fields in the data structure called
%                .lftHist   
%                .lftHist_InRegion
%                .lftHist_OutRegion
%           which contain the lifetime histogram of the objects 
%                 - in the entire image
%                 - only inside the segmented region
%                 - only outside the segmented region
% REMARKS 
%
% Dinah Loerke, last modified April 20, 2008



% number of entries in data structure
lens = length(data);

od = cd;

for i=1:lens
    
    fprintf('movie #%02d',i);
    
    % current path
    path = data(i).source;
    
    % number of frames for this exp
    lenf = data(i).movieLength;
    
    
    % load segmentation image from specified location
    SegmFileName = data(i).segmentDataFileName;
    SegmFilePath = data(i).segmentDataFilePath;

    cd(SegmFilePath);
    SegmentMask = imread(SegmFileName);    
    
    % load lifetime inof data file
    lftpath = [path,'/LifetimeInfo'];
    cd(lftpath);
    lftname = 'lftInfo.mat';
    loadfile = load(lftname);
    cd(od);
    
    lftInfo = loadfile.lftInfo;

    % calculate segmentation status (1=Inside, 0=oustide segmented region)
    [segmentStatusVector] = calcIORegionLfthistSimple(lftInfo, SegmentMask);
    
    lftMat = lftInfo.Mat_lifetime;
    statMat =  lftInfo.Mat_status;

    [sx,sy] = size(lftMat);

    lftVec = nan*zeros(sx,1);
    lftVecIN = nan*zeros(sx,1);
    lftVecOUT = nan*zeros(sx,1);

    % a trajectory is counted for the lifetime analysis if the status of 
    % the trajectory is ==1 and the value of the gaps is ==4
    for k=1:sx
        % current status vector
        cstat = nonzeros(statMat(k,:));
        % current lifetime vector
        clft = lftMat(k,:);
        if ( (min(cstat)==1) & (max(cstat)<5) )
            lftVec(k) = max(clft);
            % count into either IN or OUT vector depending on
            % segmentation status
            if segmentStatusVector(k)==1
                lftVecIN(k) = max(clft);
            else
                lftVecOUT(k) = max(clft);
            end

        end
    end    
    
    % lifetime vectors containing all counted lifetime lengths
    lftVec = lftVec(find(isfinite(lftVec)));
    lftVecIN = lftVecIN(find(isfinite(lftVecIN)));
    lftVecOUT = lftVecOUT(find(isfinite(lftVecOUT)));

    % lifetime histograms
    data(i).lftHist = hist(lftVec,[1:lenf]);
    data(i).lftHist_InRegion = hist(lftVecIN,[1:lenf]);
    data(i).lftHist_OutRegion = hist(lftVecOUT,[1:lenf]);
        
    fprintf('\b\b\b\b\b\b\b\b\b');       

end % of for

fprintf('\n'); 

end % of function
