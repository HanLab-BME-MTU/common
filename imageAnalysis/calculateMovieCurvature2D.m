function [sKappa, kappaPos, kappaNeg] = calculateMovieCurvature2D(MD,iChanInt)
% function [] = calculateMovieCurvature2D(MD,iChanInt) calcuates curvature
% of 2D mask of single-frame movie. 

% Show mask
iMask = MD.getProcessIndex('MaskRefinementProcess');
maskProc = MD.getProcess(iMask);
iChans = maskProc.checkChannelOutput;
% iChanInt = 2;
iFrame = 1;
curMask = maskProc.loadChannelOutput(iChanInt,iFrame);

% Get the vertices info
[B,~,nBD]  = bwboundaries(curMask,'noholes');
py = B{1}(:,1);
px = B{1}(:,2);
t = (1:length(px))';
t2 = [2:length(px) 1]';

% plot(px,py,'r')

% Get the curvature
% kappa = curvature(px,py,'polynom', 10);
Lines = [t t2];
Vertices = [px py];
% Need to smooth the line because the curvature estimation is from 3
% neighbors.
ppX= csaps(t, px, 0.5);
ppY= csaps(t, py, 0.5);
smX=ppval(ppX,t);
smY=ppval(ppY,t);
smVertices=[smX smY];

% kappa = LineCurvature2D(Vertices,Lines);
sKappa = LineCurvature2D(smVertices,Lines);
% curvature plot
h=figure; imshow(curMask);
hold on
scatter(smX,smY,5,sKappa, 'filled')
BlueBlackRedColorMap;
caxis([-0.3 0.3])
colorbar
savefig(h,[MD.getPath filesep 'edgeCurvatureScatter.fig'])
%Since there are both peaks and valleys, we need to report positive
%curvatures and negative ones separately.
% cmap = cmapping(kappa);
% cmapSm = cmapping(sKappa);
hp=figure;
% imshow(MD.channels_(iChanInt).loadImage(iFrame),[]), hold on
% imshow(curMask); hold on
% patch([px' nan],[py' nan],[kappa' nan],'edgecolor','interp')
% hold on
patch([smX' nan],[smY' nan],[sKappa' nan],...
    'edgecolor','interp','LineWidth',3);

BlueBlackRedColorMap;
caxis([-0.3 0.3])
colorbar
savefig(hp,[MD.getPath filesep 'edgeCurvaturePatch.fig']);

kappaPos = sKappa(sKappa>0);
kappaNeg = sKappa(sKappa<0);

disp(['Mean positive curvature: ' num2str(mean(sKappa(sKappa>0))) ' (1/px)'])
disp(['Mean negative curvature: ' num2str(mean(sKappa(sKappa<0))) ' (1/px)'])

