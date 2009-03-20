function [] = extractIntensityMat(data, dvector, tbuffer, channel, shiftvar, usemask, filenameRoot); 
% read intensities at specified positions (e.g. at the locations of
% detected and tracked objects) from image series and save to a parameter file
% 
% SYNOPSIS extractIntensityMat(data, dvector, tbuffer, channel, shiftvar, usemask, filenameRoot);  
%
% INPUT:    data:       data structure, which contains the location of the
%                       movie data in the .source field
%           dvector:    distance vector for image readout (i.e. the
%                       intensity in the image is read out up to the
%                       distance dvector from the specified object
%                       locations)
%           tbuffer:    time buffer for reading intensity before appearance
%                       and after disappearance (i.e. intensity is read
%                       e.g. 10 frames before and 10 frames after the
%                       actual detected trajectory)
%           channel (optional): if the data structure contains a field that
%                       contains the location of the directory with the 
%                       image data from which the intensities are read, 
%                       specify that field here; the function will then
%                       automatically point you to this directory and save
%                       you lots of unnecessary mouse clicking
%                       default is the field .source
%           shiftvar (optional): indicates whether or not the positions
%                       should be shifted by the vector value indicated in
%                       the field .colorShiftVector; DEFAULT = no
%                       NOTE: This functionality is relevant if you are
%                       reading intensities from a dual-channel movie,
%                       where the objects are detected in one channel and
%                       you are reading intensities or parameters from a
%                       second channel that is shifted by a few pixels
%           usemask (optional): indicate whether or not a cell outline mask
%                       should be used for the reference points, provided
%                       it's available
%                       DEFAULT = no
%           filenameRoot (optional): root of the filename for results files
%                       that are written into your source directory, e.g.
%                       'ClathrinInt' will cause a file 'ClathrinInt.mat'
%                       and a file 'ClathrinInt_Ref.mat' to be written as
%                       results files
%                       DEFAULT = 'parameterMat.mat' and 'parameterMat_Ref.mat'
    
%
% OUTPUT:   file with results is written into source directory; format of
%           the results file is a matrix with n rows and k columns, where n
%           is the number of objects/trajectories and k is the number of
%           frames and the matrix contains the intensities read from the
%           specified images at the location of the objects in the matrix
%           NOTE: this matrix has the same size as the lifetime
%           matrix, and has corresponding object rows
%
% last modified: Dinah Loerke, March 17, 2009
%


%% set reference and default values
od = cd;

% set pixel shift
shift = 0;
if nargin>4
    if shiftvar==1
        shift=1;
    end
end

% set mask use
umaskVar = 0;
if nargin>5
    if usemask==1
        umaskVar = 1;
    end
end

% set names for results files
defName = 'parameterMat.mat';
defName_ref = 'parameterMat_Ref.mat';
if nargin>6
    if ~isempty(filenameRoot)
        defName = [filenameRoot,'.mat'];
        defName_ref = [filenameRoot,'_Ref.mat'];
    end
end


%% determine image files from which intensities are read
for i=1:length(data)
    
    % select first image for the intensity channel
    if nargin>3
        ipath = getfield(data,{1,i},channel);
    else
        ipath = data(i).source;
    end
    
    % change to specified directory and look for files; files can be either
    % image files or parameter matrices
    cd(ipath);         
    [imageName, imagePath] = uigetfile({'*.tif';'*.mat'},['Select first intensity image or parameter mat in movie #',num2str(i)]); 
    
    % read complete list of images or parameter files
    completeImageName = strcat(imagePath, imageName);
    imageStackList = getFileStackNames(completeImageName);
    % and store them for later use
    ImageStackList(i).list = imageStackList;
    ImageSize(i,:) = data(i).imagesize;
    
end

% pause to allow the function to close the ui window
pause(0.1);


%% read intensities at designated positions (the positions of detected and
%% tracked objects stored in the lifetimematrix)
for i=1:length(data)
    
    % print processing update 
    fprintf('movie #%02d',i);
    fprintf('\n');
    
    % determine tracking data path, and make mpm of all positions in the
    % lifetime matrix (using a number of buffer frames before and after the
    % end of trajectories)
    trackinfopath = [data(i).source,filesep,'LifetimeInfo'];
    trackinfopath2 = data(i).source;
    MPMpos_obj = extractAllPosSideBuffer(trackinfopath, tbuffer);
    
    % determine reference positions - many intensities or parameters
    % require a reference value within the cell; for this purpose, load the 
    % reference image, if one exists, and choose random positions
    cpath = [data(i).source,filesep,'SubregionsMask'];
    if (exist(cpath)==7) & (umaskVar==1)
        cd(cpath);
        image = imread('mask0001.tif');
        mask = image';
        MPMpos_ref = extractAllPosSideBuffer_ref(MPMpos_obj, mask);
    else
        MPMpos_ref = extractAllPosSideBuffer_ref(MPMpos_obj,ImageSize(i,:));
    end
    
    % put object and reference positions together, and make time vectors
    sy = size(MPMpos_obj,2);
    MPMall_pos = [ MPMpos_obj(:,:,1) ; MPMpos_ref(:,:,1)];
    tMat_t1 = [ MPMpos_obj(:,1:2:sy,2) ; MPMpos_ref(:,1:2:sy,2) ];
    tMat_t2 = [ MPMpos_obj(:,1:2:sy,3) ; MPMpos_ref(:,1:2:sy,3) ];
    
    % if the data originate from different channels that are shifted, e.g.
    % from a red and green color channels are shifted a few pixels, then
    % that shift information for this particular movie should be stored in
    % a fiels called .colorShiftVector
    % If such a field exists, then the object positions in MPMglobal have 
    % to be shifted by the shift vector to read the appropriate position in
    % the intensity/parameter images
    % NOTE: What do the dimensions of the shiftvector mean? For example, a 
    % shift of shiftvec=[-10,-5]) means that image2 is shifted in such a 
    % way that the point im1(1,1) in image1 (clathrin) overlays point 
    % im2(11,6) in image2 (actin). Thus, the shift has to be subtracted
    % from the positions in the clathrin channel to obtain the correct
    % coordinates in the actin channel. Positions outside the shifted image
    % dimensions have to be set to nan
    
    if isfield(data,'colorShiftVector') & (shift==1)
        if ~isempty(data(i).colorShiftVector)
            [msx,msy,msz] = size(MPMall_pos);
            [isx,isy] = size(image);
            shiftx = data(i).colorShiftVector(1);
            shifty = data(i).colorShiftVector(2);
            MPMshiftx = MPMall_pos(:,1:2:msy)-shifty;
            MPMshifty = MPMall_pos(:,2:2:msy)-shiftx;
                       
            badpos = find( (MPMshiftx>isy) | (MPMshiftx<1) | (MPMshifty>isx) | (MPMshifty<1));
            MPMshiftx(badpos) = nan;
            MPMshifty(badpos) = nan;
                  
            MPMall_pos(:,1:2:msy,1) = MPMshiftx;
            MPMall_pos(:,2:2:msy,1) = MPMshifty;
        end
    end
    
    % read current images
    currImageStackList = ImageStackList(i).list;    
    total_frame_num = length(imageStackList);
    
%     for k=1:total_frame_num
%         
%         cimage = imread(imageStackList{k});
%         if k==1
%             intensityImageStack = zeros(size(cimage,1),size(cimage,2),total_frame_num);
%         end
%         intensityImageStack(:,:,k) = cimage;
%         
%     end
    
    nf = min(total_frame_num,size(MPMall_pos,2)/2);
    
    % read intensity images matrix
    iMat = extractIntensity_MPMfromImageStack(MPMall_pos, currImageStackList, dvector);
    
    % separate the imat results back into the objects and the reference
    % points
    bpoint = size(MPMpos_obj,1);
    epoint = size(MPMall_pos,1);
    
    iMat_obj = iMat(1:bpoint,1:nf,:);
    iMat_obj(:,:,2) = tMat_t1(1:bpoint,1:nf);
    iMat_obj(:,:,3) = tMat_t2(1:bpoint,1:nf);
    
    iMat_ref = iMat(bpoint+1:epoint,1:nf,:);
    iMat_ref(:,:,2) = tMat_t1(bpoint+1:epoint,1:nf);
    iMat_ref(:,:,3) = tMat_t2(bpoint+1:epoint,1:nf);    
    
    % save results into the .source directory under the specified names
    cd(trackinfopath2);
    save(defName,'iMat_obj');
    save(defName_ref,'iMat_ref');
    
    
end % of for i-loop

% return to original directory
cd(od);



end % of function



    
    