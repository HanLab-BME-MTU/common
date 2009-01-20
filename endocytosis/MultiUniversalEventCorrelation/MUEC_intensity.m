function [result] = MUEC_intensity(data, restvector, dvector, tvector, channel); 
% Multichannel Universal Event Correlator
% 
% SYNOPSIS [result] = MUEC_intensity(data, restvector, dvector, tvector); 
%
% INPUT:    data:   data structure
%           restvector: restriction vector; if restvector is a string
%           instead of a vector, then the function will try to read the
%           MPP from the structure field of the corresponding name
%           dvector: distance vector         
%           tvector: time shift vector
%
% OUTPUT:   results:        
%
% last modified: Dinah Loerke, August 28, 2008
%
% NOTE: the assumption is that the entered tvector (e.g. [-5:10]) is in
% frames and corresponds to the fastest framerate present in the analyzed 
% data (fframerate); data with slower framerates are also collected for the
% same number of frames pre- and post-event, but in the final averaging
% analysis, the data are extrapolated only for the time range of the
% fastest rate


od = cd;
nr = 10;
refplane = find(tvector==0);

if isfield(data,'framerate')
    scvec = data(1).framerate;
else
    scvec = 1;
end

for i=1:length(data)
    
    % select first image for the second (intensity) channel
    if nargin>4
        path2 = getfield(data,{1,i},channel);
    else
        path2 = data(i).channel2;
    end
    
    % change to source directory
    cd(path2);         
    [imageName, imagePath] = uigetfile('.tif',['Select first original intensity image for movie #',num2str(i)]); 
    
    completeImageName = strcat(imagePath, imageName);
    imageStackList = getFileStackNames(completeImageName);
    
    Ch2stackList(i).list = imageStackList;
    
end


for i=1:length(data)
    
    fprintf('movie #%02d',i);
    fprintf('\n');
    
    % load trackinfo file
    cpath = [data(i).source,filesep,'TrackInfoMatrices'];
    trackinfopath = [data(i).source,filesep,'LifetimeInfo'];
    cd(cpath);
    loadfile = open('trackInfo.mat');
    trackInfo  = loadfile.trackInfo;
    framerate = data(i).framerate;
    
    % load reference image, if one exists
    cpath = [data(i).source,filesep,'SubregionsMask'];
    cd(cpath);
    image = imread('mask0001.tif');
    mask = image';
       
    if isstr(restvector)
        restrictions = getfield( data(i),restvector );
    else
        % extract positions with desired restricition properties
       restrictions = restvector;
    end


    MPMglobal = CorrelateData2Pos_extractMPM2L(trackinfopath, restrictions, tvector, framerate, 1);
    
    % if the red and green color channels are shifted - i.e. if a field 
    % colorShiftVector exists - then the clathrin pit event positions in
    % MPMglobal have to be shifted by the shift vector
    % NOTE: What do the dimensions of the shiftvector mean? For example, a 
    % shift of shiftvec=[-10,-5]) means that image2 is shifted in such a 
    % way that the point im1(1,1) in image1 (clathrin) overlays point 
    % im2(11,6) in image2 (actin). Thus, the shift has to be subtracted
    % from the positions in the clathrin channel to obtain the correct
    % coordinates in the actin channel. Positions outside the shifted image
    % dimensions have to be set to nan
    
    if isfield(data,'colorShiftVector')
        if ~isempty(data(i).colorShiftVector)
            [msx,msy,msz] = size(MPMglobal);
            [isx,isy] = size(image);
            shiftx = data(i).colorShiftVector(1);
            shifty = data(i).colorShiftVector(2);
            MPMshiftx = MPMglobal(:,1:2:msy)-shifty;
            MPMshifty = MPMglobal(:,2:2:msy)-shiftx;
                       
            badpos = find( (MPMshiftx>isy) | (MPMshiftx<1) | (MPMshifty>isx) | (MPMshifty<1));
            MPMshiftx(badpos) = nan;
            MPMshifty(badpos) = nan;
                  
            MPMglobal(:,1:2:msy,1) = MPMshiftx;
            MPMglobal(:,2:2:msy,1) = MPMshifty;
        end
    end
            
            
    pause(0.1);
    
    % compact MPMglobal (get rid of nan rows)
    projpos_global = nanmean(MPMglobal(:,:,1),2);
    fpos_global = find(isfinite(projpos_global));
    MPMglobal = MPMglobal(fpos_global,:,:);
        
    % make random positions for reference measurement
    MPMbas = MPMglobal(:,:,1);
    goodpos = find( MPMglobal(:,:,2) == refplane );
    MPMbas_use = nan*MPMbas;
    MPMbas_use(goodpos) = MPMbas(goodpos);
    MPMbas = [];
    
    for k=1:nr
        MPMrandom_curr = makeRandomMPM(MPMbas_use, mask, 1);
        if k==1
            [srx,sry] = size(MPMrandom_curr);
            MPMrandom_all = zeros(10*srx,sry);
        end
        MPMrandom_all((k-1)*srx+1:k*srx,:) = MPMrandom_curr;
    end
       
    % time-shift random positions
    MPMrandomshift = CorrelateData2Pos_timeshiftMPM2L(MPMrandom_all, tvector);
    MPMrandomshift(MPMrandomshift==0) = nan;
            
    [lx1,ly1,lz1] = size(MPMglobal);
    [lx2,ly2,lz2] = size(MPMrandomshift);
    
    MPMall = [MPMglobal ; MPMrandomshift];
    
        
    imageStackList = Ch2stackList(i).list;    
    total_frame_num = length(imageStackList);
    
    for k=1:total_frame_num
        
        cimage = imread(imageStackList{k});
        if k==1
            intensityImageStack = zeros(size(cimage,1),size(cimage,2),total_frame_num);
        end
        intensityImageStack(:,:,k) = cimage;
        
    end
    
    
    CORRresults_all = CorrelateData2Pos(MPMall, intensityImageStack, dvector, 'intensity');
    
    % extract kinetic scores at desired positions
    %CORRresults_pos = CorrelateData2Pos(MPMglobal, data(i).ksnamelist, dvector, 'kscores');
    
    % extract kinetic scores for random positions
    %CORRresults_rand = CorrelateData2Pos(MPMrandomshift, data(i).ksnamelist, dvector, 'kscores');
    
    % global results
    %ResultsGlobal(:,:,i) = CORRresults_pos;
    %ResultsRandom(:,:,i) = CORRresults_rand;
    
    ResultsGlobalData(i).intmatrix = CORRresults_all(1:lx1,:,:);
    ResultsRandomData(i).intmatrix = CORRresults_all(lx1+1:lx1+lx2,:,:);
    ResultsNumTraj(i) = lx1;
        
end % of for i-loop

cd(od);


%% the raw data is now processed for averaging

% average results taking into account possible variable framerates

for i=1:length(data)
    cfr = data(i).framerate;
    if i==1
        ffr = cfr;
    else
        ffr = min(ffr,cfr);
    end
    framerates(i) = cfr;
end
% ffr ist now fastest framerate in data set


figure; hold on;

for i=1:length(data)
    
   
    ResultsGlobalCurr       = ResultsGlobalData(i).intmatrix;
    ResultsRandomCurr       = ResultsRandomData(i).intmatrix;
    ResultsRandomCurrAV     = nanmean(ResultsRandomCurr,1);
    ResultsRandomCurrAVMAT  = repmat(ResultsRandomCurrAV,ResultsNumTraj(i),1);
    ResultsGlobalCurrCORR   = ResultsGlobalCurr - ResultsRandomCurrAVMAT;
    
    % interpolate for fastest framerate is necessary
    if framerates(i)>ffr
        tff = ffr * tvector;
        tcf = framerates(i) * tvector;
        
        ResultsGlobalCurr_EP = interp1(tcf,ResultsGlobalCurr',tff);
        ResultsRandomCurrAV_EP = interp1(tcf,ResultsRandomCurrAV,tff);
        ResultsRandomCurrAVMAT  = repmat(ResultsRandomCurrAV_EP,ResultsNumTraj(i),1);
        ResultsGlobalCurrCORR   = ResultsGlobalCurr_EP' - ResultsRandomCurrAVMAT;
    end      
 
    
    
    if i==1
        ResultsGlobalAll = ResultsGlobalCurrCORR;
    else
        ResultsGlobalAll = [ResultsGlobalAll;ResultsGlobalCurrCORR];
    end
    
    [rx,ry,rz] = size(ResultsRandomCurrAV);
    ResultsGlobalCurrAV = nanmean(ResultsGlobalCurr,1);
    
    if rz>1
        for r=1:rz
            plot(scvec*tvector,ResultsRandomCurrAV(:,:,r),'r-');
            plot(scvec*tvector,ResultsGlobalCurrAV(:,:,r),'b-');
        end
    else
        plot(scvec*tvector,ResultsRandomCurrAV,'r-');
        plot(scvec*tvector,nanmean(ResultsGlobalCurr,1),'b-');
    end
    
end






% average the results of individual movies
ResultsGlobalAVE = nanmean(ResultsGlobalAll,1);
ResultsGlobalError = nanstd(ResultsGlobalAll,[],1)/sqrt(sum(ResultsNumTraj));

figure
[rx,ry,rz] = size(ResultsGlobalAVE);
if rz==1
    errorbar(scvec*tvector,ResultsGlobalAVE,ResultsGlobalError);
else
    for r=1:rz
        errorbar(scvec*tvector,ResultsGlobalAVE(:,:,r),ResultsGlobalError(:,:,r));
    end
end

amin = min(scvec*tvector)-0.5;
amax = max(scvec*tvector)+0.5;
axis([amin amax -2.5e-7 2.5e-7]);
xlabel('time point relative to CCP event (sec)');
ylabel('relative kinetic score');

% display Results
result.parfield = ResultsGlobalAVE;
result.parerror = ResultsGlobalError;

result.restrictions = restvector;
result.distvector = dvector;
result.timevector = tvector;
result.numtraj = sum(ResultsNumTraj);

end % of function



    
    