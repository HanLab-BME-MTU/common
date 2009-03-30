function [posvec] = extractPos_vecFromLftRest(lftinfopath, restvector, fr);
% this function extracts an MPM of points of interest from trackInfo using
% specified restrictions, and time delays in the third dimension
% 
% SYNOPSIS [MPMglobal] = CorrelateData2Pos_extractMPM(trackinfo,restvector, tvec, ttype);
%
% INPUT     lftinfopath       
%
%           restvector      = restriction vector, defined by
%               dstat       = restvector(1) = status [1,2,3];
%               dapp        = restvector(2) = disappearance [-1,0,1];
%               dminfr      = restvector(3) = min lifetime in frames;
%               minlft      = restvector(4) = min lifetime in seconds;;
%               maxlft      = restvector(5) = min lifetime in seconds;;
%               minlft_fr   = round(minlft/fr);
%               maxlft_fr   = round(maxlft/fr);
%
%
%           fr              = framerate; default 1
%
%
% OUTPUT:   posvec       = result MPM
%
% last modified: Dinah Loerke   01/30/2009




%% determine framerate
framerate = 1;
if nargin>2
    if ~isempty(fr)
        framerate = fr;
    end
end


%% calculate lft matrices
cd(lftinfopath);
loadfile = open('lftInfo.mat');
lmat = loadfile.lftInfo;

mat_lft   = full(lmat.Mat_lifetime);
mat_stat  = full(lmat.Mat_status);
mat_x     = full(lmat.Mat_xcoord);
mat_y     = full(lmat.Mat_ycoord);
mat_da    = full(lmat.Mat_disapp);

vec_lft   = max(mat_lft,[],2);
mat_stat(mat_stat==0) = nan;
vec_stat = nanmin(mat_stat,[],2);

% frame number
[nx,nf] = size(mat_lft);
    

%% =====================================================================
% DESIRED CONDITIONS
% for all frames in the movie, collect those locations of points that
% fulfill a number of requirements
% specifically status, minimum/maximum lifetime
[rx,ry] = size(restvector);
if min(rx,ry)>1
    nres = size(restvector,2);
else
    nres = 1;
end

for r=1:nres
    
    if nres>1
        crestvector = restvector(r,:);
    else
        crestvector = restvector;
    end
       
    % desired conditions
    dstat       = crestvector(r,1);
    dminfr      = crestvector(r,3);
    minlft      = crestvector(r,4);
    maxlft      = crestvector(r,5);
      
    minlft_fr   = round(minlft/framerate);
    maxlft_fr   = round(maxlft/framerate);
    
    %% =====================================================================

    % find positions that fulfill required conditions, as logical matrices
    % correct status
    if isfinite(dstat)
        findpos_stat = ( vec_stat==dstat );
    else
        findpos_stat = isfinite(vec_stat);
    end

    % correct lifetime    
    findpos_lft = ( (vec_lft>dminfr) & (vec_lft>=minlft_fr) & (vec_lft<=maxlft_fr) );
         
    % combine positions
    findpos_all = find(findpos_stat & findpos_lft );
    
    if nres>1
        posvec(r).posvec = findpos_all;
    else
        posvec = findpos_all;
    end
    
end % of for r-loop


end % of function

