function [msd] = MSDfromMPM_general(mpm)
% calculate matrix of mean squared displacement vectors for all
% trajectories contained in an mpm-file
% SYNOPSIS: [msd] = MSDfromMPM(mpm)
%
% INPUT:    mpm     = MPM file, where columns represent
%                     [x1 y1 x2 y2 ... xn yn]
%                     and the rows represent the trajectories
%                     NOTE that this assumes there's only one trajectory
%                     per row!!!
% OUTPUT:   msd     = mean squared displacement
%                     where every row contains the complete msd vector for
%                     every trajectory, and the columns correspond to the
%                     length of the time lag in frams
% 
% NOTE FOR POSTPROCESSING:
% NOTE 1: In this function, the msd is calculated for the maximum possible
% lag time, i.e. the length of the movie. MSD analysis is generally
% considered acceptable for lags up to HALF of the length of the actual 
% trajectory, because of decreasing statistics and increasing error for
% longer time lags. The user of the function is himself/herself responsible
% for restricting the resulting MSD length to the correct length
%
% NOTE FOR PREPROCESSING:
% NOTE 2: This function assumes that missing values in trajectories - e.g.
% as caused by gap-closing - are represented by nan values. In some
% applications, gap-closed trajectory positions are filled with the last
% available x or y-position. These would then be erroneously counted as 
% true object positions, so make sure that gap-closed positions are filled
% with NAN before running this function
%
% last modified: Dinah Loerke 03/27/2009 


[mpmx,mpmy]     = size(mpm);
frames          = round(mpmy/2);
ntau            = 1:frames-1;

    
% initialize msd matrix
msd  = nan*zeros(mpmx,length(ntau));

% extract x and y matrices from mpm
xmat = mpm(:,1:2:mpmy);
ymat = mpm(:,2:2:mpmy);

% loop over lag times
for n=1:length(ntau)
    
    fprintf('frame %04d',n);
    
    f1 = n+1;
    f2 = frames-n;
    
    crop1x = xmat(:,f1:frames);
    crop2x = xmat(:,1:f2);
    
    crop1y = ymat(:,f1:frames);
    crop2y = ymat(:,1:f2);
    
    % x- and y-displacement matrices
    cdiffx = crop1x - crop2x;
    cdiffy = crop1y - crop2y;
    
    % squared displacement matrix
    diff_squared = cdiffx.^2 + cdiffy.^2;   
    
    % average squared displacement over all available values in each traj
    msd(:,n) =  nanmean(diff_squared,2);
    
    fprintf('\b\b\b\b\b\b\b\b\b\b');
    
end % of for

end  % of function

