function [] = runDiffusionAnalysis(rest);

% runDiffusion analyses diffusion of tracks for a given condition that fit 
% the characteristics imposed by the restriction vector
%
% INPUT
%           restrict  = restriction vector can have variable length;
%                       minimum length is five, where the entries are
%                       [stat da minfr minlft maxlft]          
%
% OUTPUT
%
% REMARKS   The function does not overwrite previous results; it uses
%           secureSave. It makes a new folder called diffusionAnalysis
%           under movie folder. It uses the function
%           trackDiffusionAnalysis1 to do the diffusion analysis. The
%           restriction vector is saved alond with the diffusion analysis.
%
% Uses: trackDiffusionAnalysis1
%       loadIndividualMovies
%       plotPopulationProbabilityDistribution
%
% Daniel Nunez, July 9, 2008

%% EXPLANATION of restriction values:
% rest = [stat da minfr minlft maxlft minint maxint minmot maxmot]
%                     
% stat  =   object status, possible values are 1, 2, 3
%           1: object appears and disappears in the movie
%           2: object is present for entire length of the movie
%           3: object is present in either first or last frame of the movie
% da    =   appearance/disappearance status, possible values are 1,0,-1
%           1: marks first frame
%           0: marks middle frame
%           -1: marks last frame
% minfr =   minimum lifetime/length in FRAMES - e.g. 4, to exclude tracking
%           artifacts of false detection positives
% minlft =  minimum lifetime in SECONDS - e.g. 60 to select for productive
%           clathrin-coated pits
% maxlft =  maximum lifetime in SECONDS - e.g. 25 to select for abortive
%           clathrin-coated pits
%
% These latter criteria allow you to select e.g. the brightest 10% of the
% population, or the faster 50%.
% 
% EXAMPLE:  to select the positions where productive pits appear
%           rest1 = [1 1 4 60 300]
%           to select the positions where abortive pits are located in each
%           frame
%           rest2 = [1 0 4 8 25]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%LOAD MOVIES
[experiment] = loadIndividualMovies();

%FOR EACH MOVIE
for iexp = 1:length(experiment)

    %print movie number
    display(['running movie number ' num2str(iexp) ' out of ' num2str(length(experiment))]);

    %LOAD TRACKS
    cd(experiment(iexp).source);
    cd('TrackInfoMatrices');
    load trackInfo.mat;
    if exist('trackInfo')
        trackInfoMat = trackInfo;
    end


    %GET TRACKS OF INTEREST FROM LIFETIME ANALYSIS RESULTS
    %load lifetime analysis results
    cd([experiment(iexp).source filesep 'LifetimeInfo'])
    load('lftInfo')

    % status matrix
    statMat = full(lftInfo.Mat_status);
    % lifetime matrix
    lftMat = full(lftInfo.Mat_lifetime);
    % disapp status matrix
    daMat = (lftInfo.Mat_disapp);
    % framerate
    framerate = experiment(iexp).framerate;
    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos1 = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
    findPos2 = find((statMat==rest(2,1)) & (daMat==rest(2,2)) &...
        (lftMat>rest(2,3)) & (lftMat>round(rest(2,4)/framerate)) & (lftMat<round(rest(2,5)/framerate)));
    findPos = [findPos1;findPos2];
    [pitID,frame] = ind2sub(size(lftMat),findPos);


    %PREPARE TRACKS FOR ANALYSIS
    %zeros for position are interpreted not as missing tracks, but as point
    %on the vbertex; these are therefore made NaNs. The other non-position
    %columns left as zeros. This assumes that zeros for positions are
    %missing tracks.
    tracks = full(trackInfoMat(pitID,:));
    tracks(tracks == 0) = nan;
    tracks(:,3:8:end-5) = 0;
    tracks(:,5:8:end-3) = 0;
    tracks(:,6:8:end-2) = 0;
    tracks(:,7:8:end-1) = 0;
    tracks(:,8:8:end) = 0;

    %RUN ANALYSIS
    alphaValues = [0.05 0.1];
    checkAsym = 1;
    [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,1,2,checkAsym,alphaValues,1);

    %SAVE RESULTS
    mkdir([experiment(iexp).source filesep 'diffusionAnalysis']);
    secureSave([experiment(iexp).source filesep 'diffusionAnalysis' filesep 'diffusionAnalysisResults' datestr(now,'yyyymmdd')],'diffAnalysisRes','errFlag','alphaValues','checkAsym','rest');
    
    %clear trackInfo so that it doesn't interfere with the trackInfo of the
    %next movie (because there are two names for this variable depending on
    %the time the movie was tracked, it is possiblt for the analysis to
    %continue with the trackInfo from a previous movie..giving repeated
    %results and skiping the movie)
    clear trackInfo;
    clear trackInfoMat;
end %of for each movie

end %of function