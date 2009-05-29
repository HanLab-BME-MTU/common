function [trackClass,mssSlope,genDiffCoef,scalingPower,normDiffCoef] ...
    = trackMSSAnalysis(tracks,probDim,momentOrders,alphaMSS)
%TRACKMSSANALYSIS classifies trajectories based on their moment scaling spectrum

%SYNPOSIS [trackClass,mssSlope,genDiffCoef,scalingPower,normDiffCoef] ...
%    = trackMSSAnalysis(tracks,probDim,momentOrders,alphaMSS)
%
%INPUT  tracks      : Matrix indicating the positions and amplitudes of the
%                     tracked features. Number of rows = number of tracks,
%                     number of columns = 8*number of time points.
%                     Each row consists of
%                     [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...].
%                     NaN is used to indicate time points where the track
%                     does not exist.
%       probDim     : Problem dimensionality. Optional. Default: 2.
%       momentOrders: Orders of moments to be calculated.
%                     Optional. Default: 0 through 6.
%       alphaMSS    : Alpha-value for classification. Can take the values
%                     0.2, 0.1, 0.05 and 0.01. One can enter one value, in
%                     which case it will be used for both confined and
%                     directed, or two values, where the first will be used
%                     for confined and the second for directed.
%                     Optional. Default: 0.1 for both.
%
%OUTPUT trackClass  : # tracks x 1 vector of track classification.
%                     Values mean the following ...
%                     0 = stalled. (NOT IMPLEMENTED YET)
%                     1 = confined Brownian.
%                     2 = pure Brownian.
%                     3 = Brownian with drift (directed).
%                     NaN = not classified.
%       mssSlope    : # tracks x 1 vector of each track's slope of the line
%                     representing moment scaling power vs. moment order.
%                     NaN indicates tracks that could not be analyzed.
%       genDiffCoef : # tracks x # orders array of generalized diffusion
%                     coefficients for every moment order considered.
%                     NaN indicates tracks that could not be analyzed.
%       scalingPower: # tracks x # orders array of powers with which moment
%                     values scale with time.
%                     NaN indicates tracks that could not be analyzed.
%       normDiffCoef: # tracks x 1 vector of each track's "normal"
%                     diffusion coefficient.
%                     NaN indicates tracks that could not be analyzed.
%
%REMARKS
%(1) Algorithm is based on Ewers et al. 2005. PNAS 102: 15110-15115 and
%Ferrari et al. 2001. Physica D 154: 111-137.
%(2) Analysis assumes that there are no kinks in the moment scaling
%spectrum curve, i.e. that the motion is strongly self-similar. Weakly
%self-similar processes will generate an MSS which is piece-wise
%continuous, hence before fitting to estimate the slope the curve must be
%chopped into smaller straight-line pieces (but this is not done).
%(3)MSS slope thresholds corresponding to alphaMSS = 0.1, 0.05 and 0.01 in
%2D case are calculated using a smoothing spline fit to the corresponding
%percentiles, which were in turn derived from a sample of 35500
%simulations. For all other conditions, the threshold is calculating using
%simple piece-wise line fit to the corresponding percentiles which were
%derived from a sample of 2000 simulations.
%
%Khuloud Jaqaman, March 2008

%% input

if nargin < 1
    disp('--trackMSSAnalysis: Please input tracks to analyze.');
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(momentOrders)
    momentOrders = 0 : 6;
end
numOrders = length(momentOrders);

if nargin < 4 || isempty(alphaMSS)
    [alphaMSSConf,alphaMSSDir] = deal(0.1);
elseif length(alphaMSS) == 1
        [alphaMSSConf,alphaMSSDir] = deal(alphaMSS);
elseif length(alphaMSS) == 2
    alphaMSSConf = alphaMSS(1);
    alphaMSSDir = alphaMSS(2);
end

%get number of tracks and frames
[numTracks,numFramesMovie] = size(tracks);
numFramesMovie = numFramesMovie / 8;

%find indices of tracks that are >= 20 frames long - do not attempt
%to calculate moments for shorter tracks
criteria.lifeTime.min = 20;
indx4diff = chooseTracks(tracks,criteria);
clear criteria

%% alpha-value for classification

%determine threshold based on alpha-value and dimensionality
switch probDim
    case 1
        switch alphaMSSConf
            case 0.2 %10th percentile and 90th percentile
                mssThreshNeg = threshMSS1D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                mssThreshNeg = threshMSS1D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                mssThreshNeg = threshMSS1D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                mssThreshNeg = threshMSS1D_p01(numFramesMovie);
        end
        switch alphaMSSDir
            case 0.2 %10th percentile and 90th percentile
                [dummy,mssThreshPos] = threshMSS1D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                [dummy,mssThreshPos] = threshMSS1D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                [dummy,mssThreshPos] = threshMSS1D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                [duumy,mssThreshPos] = threshMSS1D_p01(numFramesMovie);
        end
    case 2
        switch alphaMSSConf
            case 0.2 %10th percentile and 90th percentile
                mssThreshNeg = threshMSS2D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                mssThreshNeg = threshMSS2D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                mssThreshNeg = threshMSS2D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                mssThreshNeg = threshMSS2D_p01(numFramesMovie);
        end
        switch alphaMSSDir
            case 0.2 %10th percentile and 90th percentile
                [dummy,mssThreshPos] = threshMSS2D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                [dummy,mssThreshPos] = threshMSS2D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                [dummy,mssThreshPos] = threshMSS2D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                [dummy,mssThreshPos] = threshMSS2D_p01(numFramesMovie);
        end
    case 3
        switch alphaMSSConf
            case 0.2 %10th percentile and 90th percentile
                mssThreshNeg = threshMSS3D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                mssThreshNeg = threshMSS3D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                mssThreshNeg = threshMSS3D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                mssThreshNeg = threshMSS3D_p01(numFramesMovie);
        end
        switch alphaMSSDir
            case 0.2 %10th percentile and 90th percentile
                [dummy,mssThreshPos] = threshMSS3D_p20(numFramesMovie);
            case 0.1 %5th percentile and 95th percentile
                [dummy,mssThreshPos] = threshMSS3D_p10(numFramesMovie);
            case 0.05 %2.5th percentile and 97.5th percentile
                [dummy,mssThreshPos] = threshMSS3D_p05(numFramesMovie);
            case 0.01 %0.5th percentile and 99.5th percentile
                [dummy,mssThreshPos] = threshMSS3D_p01(numFramesMovie);
        end
end

%% memory for trajectory classification

%classification means ...
%0 = stalled
%1 = confined Brownian
%2 = pure Brownian
%3 = drift/directed
%NaN = unclassified

trackClass = NaN(numTracks,1);
mssSlope = NaN(numTracks,1);
genDiffCoef = NaN(numTracks,numOrders);
scalingPower = NaN(numTracks,numOrders);
normDiffCoef = NaN(numTracks,1);

%% moments and their scaling with time

for iTrack = indx4diff'

    %get track start and end time
    trackSEL = getTrackSEL(tracks(iTrack,:));
    startTime = trackSEL(1);
    endTime = trackSEL(2);
    numTimePoints = trackSEL(3);

    %extract track's coordinates and their standard deviations
    coordinates = [tracks(iTrack,1:8:end)' tracks(iTrack,2:8:end)' tracks(iTrack,3:8:end)'];
    coordinates = coordinates(startTime:endTime,:);
    standardDevs = [tracks(iTrack,5:8:end)' tracks(iTrack,6:8:end)' tracks(iTrack,7:8:end)'];
    standardDevs = standardDevs(startTime:endTime,:);

    %define maximum time lag for moment calculation
    maxLag = min(30,floor(numTimePoints/4));

    %calculate track moments
    trackMomentsT = calcTrackMoments(coordinates,standardDevs,momentOrders,maxLag);
    trackMoments = [trackMomentsT.momentValues];
    trackMoments = trackMoments(:,1:2:end);

    %estimate the moment scaling spectrum (MSS),
    %i.e. the scaling power for all moments
    scalingPowerT = NaN(1,numOrders);
    genDiffCoefT = NaN(1,numOrders);
    for iOrder = 1 : length(momentOrders)

        %caculate ln(lag) and ln(moment)
        lnTime = log((1:maxLag)');
        lnMoment = log(trackMoments(:,iOrder));
        
        %remove any NaNs
        indxGood = find(~isnan(lnMoment));
        lnTime = lnTime(indxGood);
        lnMoment = lnMoment(indxGood);
        
        %if there are moments to fit ...
        if length(lnMoment) > 1

            %fit a straight line in the plot of lnMoment vs. lnTime
            slParam = polyfit(lnTime,lnMoment,1);

            %get scaling power and generalized diffusion coefficient
            scalingPowerT(iOrder) = slParam(1);
            genDiffCoefT(iOrder) = exp(slParam(2)) / 2 / probDim;
            
            %if this is the 2nd moment, calculate the "normal" diffusion
            %coefficient
            if momentOrders(iOrder)==2
                options = optimset('Display','off','Jacobian','on');
                lnSlope = lsqcurvefit(@strLineFun2,1,lnTime(1:min(5,...
                    length(lnTime))),lnMoment(1:min(5,length(lnMoment))),...
                    [],[],options);
                normDiffCoefT = exp(lnSlope) / 2 / probDim;                
            end
            
        end

    end

    %keep only non-NaN scaling powers
    indxGood = find(~isnan(scalingPowerT));
    momentOrders4fit = momentOrders(indxGood);
    scalingPowerT = scalingPowerT(indxGood);
    genDiffCoefT = genDiffCoefT(indxGood);
    
    %if there are non-NaN scaling powers
    if ~isempty(scalingPowerT)

        %fit a straight line to the MSS
        slParam = polyfit(momentOrders4fit,scalingPowerT,1);

        %get the slope of the line
        mssSlopeT = slParam(1);

        %classify track as ...
        %1 = confined Brownian, if MSS slope < mssThreshNeg
        %2 = pure Brownian, if mssThreshNeg <= MSS slope <= mssThreshPos
        %3 = directed, if MSS slope > mssThreshPos
        if ~isnan(mssSlopeT)
            if mssSlopeT < mssThreshNeg(numTimePoints)
                trackClass(iTrack) = 1;
            elseif mssSlopeT > mssThreshPos(numTimePoints)
                trackClass(iTrack) = 3;
            else
                trackClass(iTrack) = 2;
            end
        end

        %save additional output information
        mssSlope(iTrack) = mssSlopeT;
        genDiffCoef(iTrack,:) = genDiffCoefT;
        scalingPower(iTrack,:) = scalingPowerT;
        normDiffCoef(iTrack) = normDiffCoefT;
        
    end

end

%% subfunction 1
function [y,d] = strLineFun2(logSlope,x)

y = logSlope + x;
d = ones(size(x));


%% thresholds

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p20(nTP)

%1D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 60 100 500];
slopeM = [0.002434794446981 0.000542795021531 0.000165329365168 0];
slopeP = [-0.001541239220102 -0.00036860800289 -0.00002638219148 0];
interseptM = [0.130233710469096 0.243753675996077 0.319246807268674 0.40191148985285];
interseptP = [0.672499246818346 0.602141373785626 0.567918792644692 0.554727696904506];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p10(nTP)

%1D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 60 150 500];
slopeM = [0.00403865723972 0.00073054296121 0.000167789786604 0];
slopeP = [-0.001863373027824 -0.000207607837099 -0.000068295975124 0];
interseptM = [0.018637170190375 0.184042884115918 0.296593519037011 0.380488412339042];
interseptP = [0.727469044719516 0.628123133276013 0.607226353979748 0.5730783664177];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p05(nTP)

%1D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 200 500];
slopeM = [0.004221335524384 0.000794743062117 0.000192611669839 0];
slopeP = [-0.002478698219027 -0.000203422663306 -0.000101280672122 0];
interseptM = [-0.059561585354062 0.146033962381945 0.266460240837543 0.362766075757204];
interseptP = [0.772420366381382 0.658656588595335 0.638228190358568 0.587587854297634];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p01(nTP)

%1D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 150 500];
turnPointsP = [20 45 100 500];
slopeM = [0.00748519566052008 0.00132812830129705 0.000291477165467116 0];
slopeP = [-0.00263804313196125 -0.000652822904393821 -0.000123920022338378 0];
interseptM = [-0.27508360480994 0.0327697631512118 0.188267433525701 0.334006016259259];
interseptP = [0.839817228313045 0.750482318072511 0.697592029866967 0.635632018697778];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p20(nTP)

%2D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 60 150 500];
slopeM = [0.002365010936411 0.000327173962804 0.000137961052668 0];
slopeP = [-0.001158501020215 -0.000170757301414 -0.000033415212947 0];
interseptM = [0.227688919307209 0.329580767987585 0.367423350014689 0.436403876348735];
interseptP = [0.634413262884249 0.575148639756167 0.554547326486174 0.537839720012491];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p10(nTP)

%2D, alpha = 0.1

%NEW

mssSlopeThreshM = [...
    0.224574197458434
    0.232559328197472
    0.240357791324024
    0.246731733302628
    0.247450701750563
    0.252840787466973
    0.25825075827193
    0.263944731727673
    0.261076605760038
    0.266197697287224
    0.270178419497252
    0.274606592489297
    0.271376194203051
    0.276085785031232
    0.279378728963397
    0.283219603898803
    0.279654403572
    0.283523320918118
    0.286623354962519
    0.290307212465952
    0.286839418118361
    0.288984857293021
    0.291657181170727
    0.294696186390912
    0.291296036625286
    0.294274998957638
    0.297400403843315
    0.300093544028548
    0.296058303558528
    0.298984599873247
    0.300837413650959
    0.30355875941062
    0.300537241041706
    0.303553899696145
    0.304643574404516
    0.306858640291308
    0.304517227216722
    0.305786651507545
    0.307694416757531
    0.309570228584845
    0.307136751429856
    0.308913831979711
    0.311123413298387
    0.312495448400378
    0.309540038446651
    0.310984660138447
    0.312990691870498
    0.31484605563284
    0.312802168261889
    0.314774386519118
    0.316575956262911
    0.318306178204258
    0.31630405351131
    0.317642463153647
    0.319191262418683
    0.320771521937496
    0.318990933603294
    0.320242585948214
    0.321678022160499
    0.323025857488633
    0.321069206265912
    0.32266446190648
    0.32411447829162
    0.325132246686479
    0.322992583642119
    0.323812952639433
    0.32491153817459
    0.325757675440457
    0.323942358143971
    0.324935715890395
    0.326010916798555
    0.32679931630416
    0.324219560510537
    0.325508461984927
    0.326917158107201
    0.328018604003
    0.326232160746548
    0.326961490754226
    0.328028478958878
    0.329102484331661
    0.326973602201008
    0.327895273457125
    0.328838467242407
    0.330114641324494
    0.328381938906311
    0.329645800011516
    0.330929040311086
    0.33169733354709
    0.330143207685598
    0.331222331707293
    0.332107477065053
    0.333050657036604
    0.331447652894312
    0.332245520588352
    0.333301465387876
    0.334480220741545
    0.33275447586093
    0.333340247514732
    0.334356092454392
    0.334701601312588
    0.33305602342014
    0.333571439651457
    0.33452814107543
    0.335136860203542
    0.335675609552624
    0.336156021072486
    0.336796874380794
    0.337566147391449
    0.338197778725749
    0.338805050847811
    0.339452759437724
    0.340165806025131
    0.340815447803009
    0.34137901878133
    0.341845378660267
    0.342643968821759
    0.3430739292543
    0.343628709867637
    0.344030449547619
    0.344820345837317
    0.345385150217369
    0.346192542181834
    0.347150568143299
    0.347177550210086
    0.347851444610726
    0.348210649507165
    0.348983376323579
    0.349605264937159
    0.350069938939118
    0.350855513247267
    0.351711131132179
    0.373419649246441
    0.387597688801557
    0.39757229927659
    0.405899172670596
    0.411717331218285
    0.417404214745695
    0.421737294870718];

mssSlopeThreshP = [...
    0.653147355047145
    0.65132889440543
    0.649679483704307
    0.648255163538267
    0.638266876704396
    0.636400480426
    0.634316151103251
    0.633292570511596
    0.628000687058611
    0.627089206944701
    0.62623067199019
    0.624941557090645
    0.621290363716996
    0.620656480982382
    0.619741271556897
    0.618677674924027
    0.615191927091456
    0.614672321955848
    0.613712500298175
    0.612670790407075
    0.610480153245027
    0.609441707448592
    0.608951952696789
    0.608301812139239
    0.606731156213953
    0.606088800305914
    0.605578390971856
    0.604654925884508
    0.603241283957719
    0.602568191842692
    0.602034107891529
    0.601502966069536
    0.600150178644757
    0.599954439030026
    0.599479698077799
    0.598948995542751
    0.598335083997344
    0.597895946426736
    0.597354513864994
    0.597248773582739
    0.596221321800998
    0.596318548882293
    0.595714725703501
    0.595371138010948
    0.594542730138966
    0.59422339966226
    0.594234483216372
    0.594023321567637
    0.59319228593867
    0.593056926275941
    0.59251157375923
    0.592865374565399
    0.592463004095721
    0.592090365270891
    0.591636608943834
    0.59150254266524
    0.591293152592306
    0.590893017494081
    0.5909590461544
    0.590757060805986
    0.59008876723898
    0.589764691633642
    0.589895615061654
    0.589696629469499
    0.588985888855342
    0.588611243655994
    0.588293218254769
    0.587979469301641
    0.58789107809767
    0.587849365301885
    0.587786718740068
    0.587565596487811
    0.587333249637887
    0.587227795585252
    0.586915982025173
    0.586659657746996
    0.586621219433631
    0.58655825918332
    0.586607820875872
    0.586607465763016
    0.58662852895352
    0.586521269181149
    0.586348328329744
    0.585993194401089
    0.58618479909948
    0.585884464694563
    0.585452261779234
    0.585354836543565
    0.585008661395476
    0.58491999976014
    0.584798659054402
    0.584778811701729
    0.584626919852783
    0.584563759450215
    0.584355849930455
    0.584438181230103
    0.584170392306253
    0.583949214171273
    0.583898991132972
    0.583504243272866
    0.583236345638302
    0.583088957376113
    0.583173276719509
    0.583381836019792
    0.583238013229126
    0.583076490844526
    0.582993336593433
    0.582836284243519
    0.582883399895619
    0.582625669153702
    0.582469259862584
    0.582276169598837
    0.581972357215994
    0.581853795824746
    0.581640888158356
    0.581613113203099
    0.581345099050163
    0.581257980303522
    0.581101202181861
    0.581106338839736
    0.581016085872333
    0.580823745913101
    0.580610008520541
    0.580534960994928
    0.580391660548403
    0.58031338183263
    0.580133695745293
    0.579942962575481
    0.579831508049998
    0.579685855341995
    0.579624446685547
    0.573956870902088
    0.569736067031574
    0.566437014156808
    0.563380261988323
    0.560726567291193
    0.558756111946602
    0.556881403376681];

%fit smoothing spline to threshold for confined, and interpolate
splineFitM = csaps([20:150 200:50:500],mssSlopeThreshM,0.05);
mssThreshNeg = [NaN(19,1); fnval(splineFitM,(20:min(500,nTP))')];
mssThreshNeg = [mssThreshNeg; mssThreshNeg(end)*ones(max(0,nTP-500),1)];

%fit smoothing spline to threshold for directed, and interpolate
splineFitP = csaps([20:150 200:50:500],mssSlopeThreshP,0.05);
mssThreshPos = [NaN(19,1); fnval(splineFitP,(20:min(500,nTP))')];
mssThreshPos = [mssThreshPos; mssThreshPos(end)*ones(max(0,nTP-500),1)];

%OLD
% %threshold curve parameters
% turnPointsM = [20 50 200 500];
% turnPointsP = [20 50 100 500];
% slopeM = [0.002834926966999 0.00040315841383 0.000171843106892 0];
% slopeP = [-0.001737592454918 -0.000337743432604 -0.000073809631405 0];
% interseptM = [0.169732430047339 0.291320857705762 0.337583919093415 0.423505472539345];
% interseptP = [0.684919949284652 0.614927498168961 0.588534118049036 0.551629302346427];
%
% %threshold curve evaluation
% [mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
%     slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p05(nTP)

%2D, alpha = 0.05

%NEW

mssSlopeThreshM = [...
    0.17851880381419
    0.1871730025546
    0.196806677145191
    0.204158488839211
    0.206780204232365
    0.21353652330535
    0.218577893947974
    0.223689904246123
    0.222886969243971
    0.228177104641715
    0.233663785913629
    0.239823425476981
    0.237133768841959
    0.241532853194731
    0.245839766102515
    0.24957179742297
    0.247191406774247
    0.25063478522449
    0.255473793871516
    0.25942938968348
    0.254931053813695
    0.257698407885778
    0.259785064837744
    0.263067478299164
    0.260139812933828
    0.263102433740553
    0.26741770863631
    0.269690336003608
    0.26690500294172
    0.269279550337408
    0.272191467568818
    0.274732878983067
    0.271025801318721
    0.273790820577884
    0.275009039118091
    0.27733986831005
    0.273170815788627
    0.27552019844457
    0.277916692959511
    0.280252910270569
    0.276984038010674
    0.280015495232136
    0.281743742377518
    0.284209925547003
    0.281004644730678
    0.283080116510187
    0.28478778026557
    0.286788305337569
    0.284136898032984
    0.285859418906338
    0.288250931610169
    0.289726410517383
    0.287774661557119
    0.290110792934316
    0.291577670301109
    0.29357715164282
    0.291152695818162
    0.292723741229145
    0.294428705734903
    0.296284348044087
    0.294079793467697
    0.296140523245772
    0.29674929600504
    0.298074416347853
    0.296021373583597
    0.297405492500614
    0.298702824136024
    0.299609274920867
    0.29748524233694
    0.298922935244787
    0.300396681779992
    0.301695357655895
    0.298810183668257
    0.29963335929538
    0.300726622866354
    0.302011538115881
    0.299680870615641
    0.300555178433349
    0.301500035531218
    0.302814475423732
    0.300495033226206
    0.301043161036347
    0.302599629593006
    0.304027996171868
    0.30248430433422
    0.30444287588424
    0.305919003051943
    0.307468262932956
    0.305461650063214
    0.305971272805794
    0.307098804618111
    0.308407242943892
    0.306061329632763
    0.307536501389644
    0.308452662919498
    0.309174325167614
    0.307522529764362
    0.30844934724475
    0.309982854436816
    0.310823818135481
    0.308077360907993
    0.308742662619394
    0.309647805603832
    0.310535528309306
    0.311521329842896
    0.312631583310933
    0.313612564360577
    0.31389812255908
    0.314777820622545
    0.315522433002325
    0.316141244189818
    0.316457655785337
    0.317268181024399
    0.31847940657704
    0.319288068829374
    0.319882469425602
    0.320512345395368
    0.321178015943994
    0.321707707451331
    0.322553495831334
    0.323082097331497
    0.323683265040748
    0.324461676032448
    0.32534199866625
    0.325718062141354
    0.326033897565515
    0.326711715628724
    0.327333834966668
    0.327771217677864
    0.328143860117064
    0.328532260659434
    0.354657628368557
    0.370827311917658
    0.382110290956925
    0.391079023933105
    0.397639366798427
    0.403733669063502
    0.40847212189647];

mssSlopeThreshP = [...
    0.687515318409466
    0.685276840381168
    0.682310561481955
    0.680269322151898
    0.670761438237813
    0.668850618151345
    0.667108714557779
    0.665257873482051
    0.659544959778331
    0.658564624033264
    0.656390823739982
    0.655334669652393
    0.649768630980414
    0.648554688106195
    0.64733206058663
    0.646667472432749
    0.643087317951418
    0.64260041945609
    0.642015155100117
    0.641235202334666
    0.638151435624193
    0.637565208573271
    0.636877368331779
    0.636296240775795
    0.633417881281766
    0.632710957101116
    0.632570137717717
    0.631709198598634
    0.630150270258322
    0.629673211058848
    0.629091453066231
    0.628097808090468
    0.626119162802377
    0.625861483439657
    0.624894318620891
    0.624988101259194
    0.62436785722323
    0.623850770564247
    0.623149817439354
    0.623059336956896
    0.622380571473926
    0.621447945819085
    0.620770711134517
    0.620261382175091
    0.618850229073762
    0.618680695343894
    0.618461413581746
    0.618349907588212
    0.617869911300813
    0.617768116997836
    0.617394090514407
    0.617088958117533
    0.617024690153813
    0.616723442005675
    0.616541996342523
    0.61591454854244
    0.616446171462252
    0.615821903302535
    0.615583614420069
    0.614911353866245
    0.615412525977224
    0.61498245966794
    0.614868983642384
    0.614661432898643
    0.614175193590041
    0.613519810762257
    0.612769342190297
    0.61184095238639
    0.611744329927023
    0.611383112637346
    0.611363426353595
    0.611205148201553
    0.610899733546305
    0.610600033612161
    0.610170428732971
    0.609781040806253
    0.609950334356465
    0.609734506809841
    0.609660766144285
    0.609319162924641
    0.609017707927956
    0.608680454051921
    0.60874235141495
    0.608744169581691
    0.608721091140155
    0.608754744897785
    0.608469360544911
    0.608222630220981
    0.608118605790019
    0.607514470023567
    0.607487244672096
    0.607103519035866
    0.606843905993827
    0.606598992613698
    0.606337998851331
    0.606284692420762
    0.606101286455309
    0.605932723503751
    0.605915261865577
    0.605487021856316
    0.605295643449468
    0.605061653341967
    0.605000933035169
    0.604544765962621
    0.604574912170413
    0.604567842780715
    0.604519690341924
    0.604670488535696
    0.604365550783769
    0.60423921775345
    0.604211047055473
    0.603968717484945
    0.60376937750576
    0.603330777601416
    0.603337530254406
    0.602819862688282
    0.602798789189528
    0.602620510303525
    0.60232890161527
    0.602155822077865
    0.601928969953038
    0.602114550246144
    0.601926971831044
    0.601756943301308
    0.601712245628021
    0.601572471487258
    0.601174723109729
    0.601094172824707
    0.600931866419486
    0.600479838986798
    0.600445339732039
    0.593025847204398
    0.587773525409125
    0.582762155913914
    0.578935961761433
    0.576435335811663
    0.573604611211088
    0.570678214442582];

%fit smoothing spline to threshold for confined, and interpolate
splineFitM = csaps([20:150 200:50:500],mssSlopeThreshM,0.05);
mssThreshNeg = [NaN(19,1); fnval(splineFitM,(20:min(500,nTP))')];
mssThreshNeg = [mssThreshNeg; mssThreshNeg(end)*ones(max(0,nTP-500),1)];

%fit smoothing spline to threshold for directed, and interpolate
splineFitP = csaps([20:150 200:50:500],mssSlopeThreshP,0.05);
mssThreshPos = [NaN(19,1); fnval(splineFitP,(20:min(500,nTP))')];
mssThreshPos = [mssThreshPos; mssThreshPos(end)*ones(max(0,nTP-500),1)];

%OLD
% %threshold curve parameters
% turnPointsM = [20 60 200 500];
% turnPointsP = [20 50 100 500];
% slopeM = [0.002405756738527 0.000449948051639 0.000201793444222 0];
% slopeP = [-0.001937063703976 -0.00055853262239 -0.000095409685751 0];
% interseptM = [0.14034289465302 0.257691415866306 0.307322337349635 0.40821905946088];
% interseptP = [0.724427540398996 0.6555009863197 0.609188692655777 0.561483849780384];
% 
% %threshold curve evaluation
% [mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
%     slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p01(nTP)

%2D, alpha = 0.01

%NEW

mssSlopeThreshM = [...
    0.0869293372917652
    0.10150310109063
    0.110458946420983
    0.120032758710204
    0.126841231786467
    0.132788690031883
    0.14111701308745
    0.148781687505228
    0.147252880610777
    0.1553171871407
    0.161442969632278
    0.168518651326214
    0.166927249661439
    0.169713556194574
    0.176200144207776
    0.183683423129845
    0.180635673565965
    0.184448961101505
    0.190868182290257
    0.195653039305515
    0.191252194839295
    0.196634636222973
    0.199200620916018
    0.204966948598132
    0.199034977778829
    0.203794068441596
    0.20708633211622
    0.210645344880564
    0.206685150305956
    0.210066380391817
    0.212767942989253
    0.217133344918783
    0.214240086206739
    0.214772194218289
    0.219544508270378
    0.223154571256961
    0.219011271077927
    0.220986503689222
    0.223488753108885
    0.22619483084937
    0.224194973560089
    0.226178225803561
    0.227406761136877
    0.22940624261505
    0.22740780991765
    0.228409204872548
    0.231202813362898
    0.233544129082427
    0.230112799676474
    0.231477087623101
    0.232771262427577
    0.235423552060439
    0.233230474179764
    0.235004551394325
    0.236837833499403
    0.239925178524418
    0.23643883383143
    0.238525271507289
    0.239494673183166
    0.241644865269114
    0.238530040907104
    0.239074346365199
    0.240429497038729
    0.24249052420283
    0.239960108432541
    0.242624556563014
    0.246548500562546
    0.2482388052189
    0.24607834810935
    0.24631986902394
    0.247705918615776
    0.248557756204567
    0.24508727005309
    0.246167649594781
    0.247958911621948
    0.250106631141259
    0.246563259798458
    0.24787643982634
    0.249487947306661
    0.251642114429472
    0.250353774128162
    0.253105520989134
    0.255154736842992
    0.255458715941156
    0.254481294256244
    0.255204010387031
    0.256705126931167
    0.257862965000293
    0.255465779313494
    0.256394450063663
    0.257255518773257
    0.257428789612193
    0.255947747739563
    0.257097813890015
    0.258153306063873
    0.259761729031251
    0.25848857662748
    0.259757136368261
    0.260472217040183
    0.261636814313246
    0.260411740732779
    0.261597017320547
    0.262143363765362
    0.263296152893485
    0.263778320432011
    0.264210048669862
    0.264881406247162
    0.265848493156194
    0.266330931199834
    0.266657906952071
    0.267913077915911
    0.26820345160375
    0.269306460578444
    0.270893072017012
    0.272270802899136
    0.272963520827481
    0.273538118818612
    0.2747511616727
    0.275010209388528
    0.276188398325692
    0.277080112179737
    0.277732438576256
    0.27797244279239
    0.279312295496952
    0.2803541422979
    0.280282682708896
    0.281410512555297
    0.282017075629052
    0.283308619337792
    0.284034900517493
    0.285391503662665
    0.316483666592547
    0.337051580502542
    0.351205409454649
    0.362210038176644
    0.370424610260052
    0.378363526295147
    0.38448542532347];

mssSlopeThreshP = [...
    0.749099569987682
    0.7473976079156
    0.743637214208967
    0.740737922870941
    0.731319375604135
    0.727351468795612
    0.725821086646462
    0.724427338813107
    0.719386888932021
    0.717229727042768
    0.715003842279892
    0.713661239816129
    0.70805040499508
    0.707575173780406
    0.704508899778484
    0.702332389387305
    0.699574001331452
    0.699606291151099
    0.698068508021262
    0.696006997617741
    0.69318723992507
    0.692211025007642
    0.690722524577197
    0.689626925037122
    0.685502965263368
    0.684181346563977
    0.684085797605733
    0.682477915849491
    0.679973394790633
    0.679208346225125
    0.678649097784419
    0.677296246097386
    0.675792489830621
    0.674689295772486
    0.673758471165635
    0.673234947975078
    0.67193116470555
    0.670740252773019
    0.669718527533529
    0.66955461847328
    0.670152182453987
    0.668977915409172
    0.668319144494854
    0.668084298910368
    0.667357093497072
    0.666388254370274
    0.66545390766008
    0.664901381763898
    0.664473841187214
    0.66383380074169
    0.663806781930918
    0.662659134136517
    0.662663410366108
    0.661975143813824
    0.66121695408175
    0.660557483414595
    0.658701307455414
    0.658318530191394
    0.657983579766347
    0.657075156133082
    0.657022613723566
    0.656310351440114
    0.65618928817467
    0.656120193055426
    0.655970827393449
    0.655720225002854
    0.655341431833804
    0.654715058668071
    0.65374795119537
    0.652704467781272
    0.652328414113147
    0.651964010931146
    0.652102231331959
    0.651595338799864
    0.65138981440537
    0.651197418745217
    0.650417323616883
    0.650134161247124
    0.649449085923729
    0.649007030471857
    0.64851028166696
    0.648626495681014
    0.648281264406385
    0.64784553574051
    0.647145897739493
    0.646541716627282
    0.646485311136555
    0.646116242821253
    0.646179150612063
    0.646008189636678
    0.645919205713778
    0.646035469268067
    0.646143872962349
    0.645460631374134
    0.645272904691685
    0.645075912479232
    0.645907654922474
    0.645137363484374
    0.645459383449963
    0.64527969781766
    0.646256769526858
    0.646695608933769
    0.646483957551041
    0.64566752089819
    0.645367142013856
    0.645195097085672
    0.645278667798355
    0.645587562820339
    0.645445730695456
    0.645223775994562
    0.645147789749733
    0.644474033445713
    0.644398112861634
    0.644541802659769
    0.64424830538345
    0.643874996092825
    0.643572221938234
    0.642951955763577
    0.642535125193078
    0.642109381986806
    0.641787294207641
    0.641336662383792
    0.640832914882883
    0.640613536204056
    0.640506272345644
    0.640334611469348
    0.640560495345698
    0.640354449811109
    0.640160366304582
    0.640200028639871
    0.639742976987077
    0.628255004260535
    0.622209478136329
    0.616361655002305
    0.611428918785076
    0.60592195201232
    0.601327582857087
    0.598222885543451];

%fit smoothing spline to threshold for confined, and interpolate
splineFitM = csaps([20:150 200:50:500],mssSlopeThreshM,0.05);
mssThreshNeg = [NaN(19,1); fnval(splineFitM,(20:min(500,nTP))')];
mssThreshNeg = [mssThreshNeg; mssThreshNeg(end)*ones(max(0,nTP-500),1)];

%fit smoothing spline to threshold for directed, and interpolate
splineFitP = csaps([20:150 200:50:500],mssSlopeThreshP,0.05);
mssThreshPos = [NaN(19,1); fnval(splineFitP,(20:min(500,nTP))')];
mssThreshPos = [mssThreshPos; mssThreshPos(end)*ones(max(0,nTP-500),1)];

%OLD
% %threshold curve parameters
% turnPointsM = [20 50 150 500];
% turnPointsP = [20 60 200 500];
% slopeM = [0.00435441385708182 0.000812259003939792 0.000273699080289468 0];
% slopeP = [-0.00172943051183377 -0.000284948443095389 -0.000153857046768591 0];
% interseptM = [-0.00562462883719271 0.171483113819909 0.252267102367457 0.389116642512191];
% interseptP = [0.779260212544801 0.692591288420498 0.666373009155138 0.589444485770843];
%
% %threshold curve evaluation
% [mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
%     slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p20(nTP)

%3D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 60 100 500];
slopeM = [0.001359327197525 0.000280732353905 0.000108640351462 0];
slopeP = [-0.000746906632969 -0.000348357342053 -0.000037074457479 0];
interseptM = [0.299405762325546 0.364121452942741 0.39853985343145 0.452860029162242];
interseptP = [0.609183824504497 0.585270867049586 0.554142578592106 0.535605349852791];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p10(nTP)

%3D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 100 500];
slopeM = [0.001587591010958 0.000334860802492 0.000129492806027 0];
slopeP = [-0.001033885285809 -0.000409579718459 -0.000060967307118 0];
interseptM = [0.258642983727984 0.333806796235989 0.374880395529001 0.439626798542322];
interseptP = [0.644618927677761 0.613403649310276 0.57854240817611 0.548058754617302];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p05(nTP)

%3D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 80 500];
slopeM = [0.001897773535643 0.000390310143108 0.00014504631569 0];
slopeP = [-0.001281266825658 -0.000815336722161 -0.000079326637385 0];
interseptM = [0.217433363211081 0.307881166763189 0.356933932246803 0.429457090091877];
interseptP = [0.683369058696754 0.660072553521921 0.601191746739831 0.56152842804728];

% threshold curve evaluation

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p01(nTP)

%3D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 90 200 500];
slopeM = [0.00330864748993062 0.000607571513122873 0.000183704012624133 0];
slopeP = [-0.0011001388082068 -0.000401685416953439 -7.31736162754304e-05 0];
interseptM = [0.100955106735165 0.236008905575552 0.3207824056753 0.412634411987367];
interseptP = [0.747201850658008 0.684341045445206 0.618638685309604 0.582051877171889];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

%% threshold curve evaluation subfunction

function [mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,...
    turnPointsP,slopeM,slopeP,interseptM,interseptP,nTP)

%confined diffusion threshold
mssThreshNeg = NaN(1,turnPointsM(1)-1);
for i = 1 : length(turnPointsM)-1
    x = turnPointsM(i) : turnPointsM(i+1)-1;
    mssThreshNeg = [mssThreshNeg slopeM(i)*x+interseptM(i)]; %#ok<AGROW>
end
x = turnPointsM(end) : nTP;
mssThreshNeg = [mssThreshNeg slopeM(end)*x+interseptM(end)];

%directed diffusion threshold
mssThreshPos = NaN(1,turnPointsP(1)-1);
for i = 1 : length(turnPointsP)-1
    x = turnPointsP(i) : turnPointsP(i+1)-1;
    mssThreshPos = [mssThreshPos slopeP(i)*x+interseptP(i)]; %#ok<AGROW>
end
x = turnPointsP(end) : nTP;
mssThreshPos = [mssThreshPos slopeP(end)*x+interseptP(end)];

