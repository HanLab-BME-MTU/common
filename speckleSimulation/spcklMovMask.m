function spcklMovMask
% apply mask to tiff series and save mask for each image as 'cell_mask' in
% 'edge' folder.  also save masked images in subdirectory of 'tiffs' folder

projDir=uigetdir(pwd,'Choose project directory');
inputDir=uigetdir(projDir,'Choose tiff image directory on which to apply mask');
outputDir=[inputDir filesep 'masked_tifs'];
if ~isdir(outputDir)
    mkdir(outputDir);
end

maskDir=[projDir filesep 'analysis' filesep 'edge' filesep 'cell_mask'];


[listOfFiles,tokenList] = searchFiles('.tif',[],inputDir);
nFrames=length(listOfFiles);
s=length(num2str(nFrames));
strg=sprintf('%%.%dd',s);

[maskL maskW]=size(double(imread([inputDir filesep listOfFiles{1,1}])));

mask=zeros(maskL,maskW);
mask(36:220,36:220)=1; 


% [X Y]=meshgrid(1:1:maskW,1:1:maskL);
% pts=[Y(:) X(:)];
% [dist]=dist2Pt(pts,[0.5*maskL 0.4*maskW]);
% mask(dist<(0.5*min([maskL maskW])-20))=1;
% mask(:,1:0.2*maskW)=0;


for i=1:nFrames
    indxStr=sprintf(strg,i);
    im=[inputDir filesep listOfFiles{i,1}];
    % mask=circshift(mask,[0 1]); % this moves the mask by a pixel each frm
    % im1=flipdim(double(imread(im)),2); % changes direction of flow
    im1=uint16(double(imread(im)).*mask);
    imwrite(im1,[outputDir filesep listOfFiles{i,1}]);
    imwrite(logical(mask), [maskDir filesep 'mask_image' indxStr '.tif']);
end


