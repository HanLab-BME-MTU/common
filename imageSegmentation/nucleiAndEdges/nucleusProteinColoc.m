function [nameSlide,numDetectedNuc,numAllNuc,ratioDetected,meanIntenCD45,meanIntenCD45ND]=nucleusProteinColoc(pathProject,pathImgDAPI,readCD45)
%[nameSlide,numDetectedNuc,numAllNuc,ratioDetected,meanIntenCD45]=nucleusProteinColoc(pathProject,pathImgDAPI,readCD45)
% detects and segments nucleii from the first
% image, and per each identified nucleus, finds if it is overlapping with
% the signal from the second image. The signal is determined by the maximum
% of the image intensity in the control image. Output is the ratio of the
% number of nucleii that overlap with protein signal to the all identified
% nucleii.
% input:
%       pathProject              : path to project
%       pathImgDAPI            : path to DAPI image
%       readCD45                : true if you have CD45 channel (default:
%                                       false)
% output:
%       nameSlide               : the name of microscopy image sample
%       numDetectedNuc     : the number of nucleii colocalized with images
%                                       in channel 3
%       numAllNuc               : the number of all segmented nucleii
%       ratioDetected          : numDetectedNuc/numAllNuc
%       meanIntenCD45       : mean intensity of channel 2(CD45) from
%                                       detected nucleii
% 
% also look at: batchNecleiiProteinColoc
% Sangyoon Han May 2015
if nargin<3
    readCD45=false;
end
%% Nucleii detection and segmentation
% import DAPI image - later I have to use uigetdir 
% pathProject = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Claudia/NucleusProteinColocalization/GB15-13A/Tra-CD45';
[patientName,traOrBstr]=fileparts(pathProject);
[~,patientName]=fileparts(patientName);
[~,nameDAPI]=fileparts(pathImgDAPI);
nameSlide{1}= strrep(nameDAPI, 'c1', '');

pathOut = [pathProject filesep nameSlide{1} '-output'];
if ~exist(pathOut','dir')
    mkdir(pathOut)
end
% pathImgDAPI = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Claudia/NucleusProteinColocalization/GB15-13A/Tra-CD45/GB15-13A-Tra-CD45008c1.tif';
imgDAPI = imread(pathImgDAPI);
% imgDAPI = imread('/project/cellbiology/gdanuser/adhesion/Sangyoon/Claudia/NucleusProteinColocalization/GB15-13A/Tra-CD45/GB15-13A-Tra-CD45008c1.tif');
% imgDAPI = imread('/project/cellbiology/gdanuser/adhesion/Sangyoon/Claudia/NucleusProteinColocalization/GB15-11A/Tra-CD45/GB15-11A-Tra-CD45c1.tif');
% figure, imshow(imgDAPI,[]), hold on
% pick the centers of each nucleus
% sigma = getGaussianPSFsigmaFromData(double(imgDAPI));
% [pstruct,mask] = pointSourceDetection(imgDAPI,sigma,'FitMixtures',true,'alpha',0.01);
% plot(pstruct.x,pstruct.y,'go')
% figure, imshow(mask)
% figure, histogram(pstruct.A,50)
pathLabelDAPI = [pathOut filesep patientName '-' traOrBstr 'LabelDAPI.tif'];
pathLabelDAPIwithDet = [pathOut filesep patientName '-' traOrBstr 'LabelDAPIwithDetection.tif'];
pathBwTra = [pathOut filesep patientName '-' traOrBstr 'traThresholded.tif'];
pathResults = [pathOut filesep patientName '-' traOrBstr 'results.csv'];
pathResultsCD45 = [pathOut filesep patientName '-' traOrBstr 'CD45intensity.csv'];
pathCD45withDet = [pathOut filesep patientName '-' traOrBstr 'CD45withDetection.tif'];
%% Run segmentNucleii
[labelDAPI,nLabel] = segmentNucleii(imgDAPI);
%% Showing splitted result:
% oldNLable = nLabel;
% nLabel=newLabel;
bg = [0 0 0];
n_colors = nLabel-1;
colors = distinguishable_colors(n_colors,bg);
colors = [bg; colors];

rgbDAPI = ind2rgb(labelDAPI,colors);
imwrite(rgbDAPI,pathLabelDAPI)
% figure, imshow(labelDAPIws,colors)
% i have to store this figure automatically... I'll work on the
% input/output system later on ...
%% Detect signal from protein channel
% import protein image of interest
if readCD45
    pathImgSignal = strrep(pathImgDAPI,'c1.tif','c3.tif');
else
    pathImgSignal = strrep(pathImgDAPI,'c1.tif','c2.tif');
end
imgSignal = imread(pathImgSignal);
% figure, imshow(imgSignal,[])
% Get the threshold value from control image
% pathImgSignalCtrl='/project/cellbiology/gdanuser/adhesion/Sangyoon/Claudia/NucleusProteinColocalization/GB15-13A/Tra-CD45/GB15-13A-Tra-CD45-ctrc3.tif';
% imgSignalControl = imread(pathImgSignalCtrl);
% figure, imshow(imgSignalControl,[])
% maxIntenFromCtrl = max(imgSignalControl(:));
% maxIntenFromSignal = max(imgSignal(:));
% figure, histogram(imgSignal(:),100)
% figure, histogram(imgSignalControl(:),100)
thresSignal = 3000;%quantile(imgSignalControl(:),0.999999); % this gives around 2500
bwSignal = imgSignal>thresSignal;
% figure, imshow(bwSignal)
% sum(bwSignal(:)) % there are about 182 of positive pixels
% Now Colocalization
detectedNucleii = labelDAPI(bwSignal);
[ySignal, xSignal] = ind2sub(size(bwSignal),find(bwSignal(:)));
% combImg(:,:,3) = labelDAPIrosin>0;
% combImg(:,:,1) = bwSignal;
% figure, imshow(double(combImg),[])
% [detectedNucleii,ia,ic]= unique(detectedNucleii);
[detectedNucleii,ia] = setdiff(detectedNucleii,0);
% figure, imshow(labelDAPI>0)
h=figure; imshow(labelDAPI,colors)
hold on, plot(xSignal(ia),ySignal(ia),'yo')
print('-dtiff', '-r300', pathLabelDAPIwithDet)
imwrite(bwSignal,pathBwTra)
close(h)
%% key statistics: number of detected nucleii vs. all nucleii
numDetectedNuc = length(detectedNucleii);
numAllNuc =nLabel;
ratioDetected = numDetectedNuc/numAllNuc;
%% Record CD45 mean intensity in detected nucleii
if readCD45
    meanIntenCD45 = zeros(numDetectedNuc,1);
    meanIntenTra = zeros(numDetectedNuc,1);
    pathCD45 = strrep(pathImgDAPI,'c1.tif','c2.tif');
    imgCD45 = imread(pathCD45);
    notDetectedNucleii = setdiff((1:numAllNuc)',detectedNucleii);
    numNotDetectedNuc = length(notDetectedNucleii);
    meanIntenCD45ND = zeros(numNotDetectedNuc,1); %ND: NotDetected
    meanIntenTraND = zeros(numNotDetectedNuc,1); %ND: NotDetected
    % Get the pixel list of detected adhesions
    regionDAPI = regionprops(labelDAPI,'PixelID');
    for ii=1:numDetectedNuc
        curPixelIDs = regionDAPI(detectedNucleii(ii)).PixelIdxList;
        curMeanIntenCD45 = mean(imgCD45(curPixelIDs));
        curMeanIntenTra = mean(imgSignal(curPixelIDs));
        meanIntenCD45(ii)=curMeanIntenCD45;
        meanIntenTra(ii)=curMeanIntenTra;
    end
    for ii=1:numNotDetectedNuc
        curPixelIDs = regionDAPI(notDetectedNucleii(ii)).PixelIdxList;
        curMeanIntenCD45 = mean(imgCD45(curPixelIDs));
        curMeanIntenTra = mean(imgSignal(curPixelIDs));
        meanIntenCD45ND(ii)=curMeanIntenCD45;
        meanIntenTraND(ii)=curMeanIntenTra;
    end
    
    meanMeanIntenCD45 = mean(meanIntenCD45);
    meanMeanIntenTra = mean(meanIntenTra);
    h2=figure; imshow(imgCD45,[])
    hold on, plot(xSignal(ia),ySignal(ia),'yo')
    print('-dtiff', '-r300', pathCD45withDet)
    close(h2)
else
    meanIntenCD45 = NaN;
    meanMeanIntenCD45 = NaN;
    meanIntenTra = NaN;
    meanMeanIntenTra = NaN;
    notDetectedNucleii=NaN;
    meanIntenCD45ND=NaN;
    meanIntenTraND=NaN;
end
%% saving
A = table(nameSlide,numDetectedNuc,numAllNuc,ratioDetected,meanMeanIntenCD45,meanMeanIntenTra);
writetable(A,pathResults)
if readCD45
    if numNotDetectedNuc>numDetectedNuc
        detectedNucleiiNaN = NaN(numNotDetectedNuc,1);
        meanIntenCD45NaN = NaN(numNotDetectedNuc,1); 
        meanIntenTraNaN = NaN(numNotDetectedNuc,1); 
        detectedNucleiiNaN(1:numDetectedNuc) = detectedNucleii;
        meanIntenCD45NaN(1:numDetectedNuc) = meanIntenCD45;
        meanIntenTraNaN(1:numDetectedNuc) = meanIntenTra;
        B = table(detectedNucleiiNaN,meanIntenCD45NaN,meanIntenTraNaN,notDetectedNucleii,meanIntenCD45ND,meanIntenTraND);
        writetable(B,pathResultsCD45)
    else
        notDetectedNucleiiNaN = NaN(numDetectedNuc,1);
        meanIntenCD45NDNaN = NaN(numDetectedNuc,1); 
        meanIntenTraNDNaN = NaN(numDetectedNuc,1); 
        notDetectedNucleiiNaN(1:numNotDetectedNuc) = notDetectedNucleii;
        meanIntenCD45NDNaN(1:numNotDetectedNuc) = meanIntenCD45ND;
        meanIntenTraNDNaN(1:numNotDetectedNuc) = meanIntenTraND;
        B = table(detectedNucleii,meanIntenCD45,meanIntenTra,notDetectedNucleiiNaN,meanIntenCD45NDNaN,meanIntenTraNDNaN);
        writetable(B,pathResultsCD45)
    end
end
