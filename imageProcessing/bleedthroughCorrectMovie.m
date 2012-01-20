function movieData = bleedthroughCorrectMovie(movieData,paramsIn)
%BLEEDTHROUGHCORRECTMOVIE corrects for bleedthrough of other fluorophores in the input movie.
% 
% movieData = bleedthroughCorrectMovie(movieData)
% 
% movieData = bleedthroughCorrectMovie(movieData,paramsIn)
% 
%
% This function corrects for bleedthrough of other fluorophores into a
% particular channel in the input movie. This is done using bleedthrough
% coefficients calculated from a "bleedthrough movie" where only one
% fluorophore is present, using processBleedthroughMovie.m. These
% coefficients are then used, in combination with images of the fluorophore
% channels which are bleeding into the image to be corrected, to remove the
% effects of bleedthrough.
% 
% 
%
% Input:
% 
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the corrected images to.
%       Corrected images for different channels will be saved as
%       sub-directories of this directory. If not input, the corrected
%       images will be saved to the same directory as the movieData, in a
%       sub-directory called "bleedthrough_corrected_images"
%
%       ('ChannelIndex'-> Positive integer scalar)
%       The integer index of the channel to perform bleedthrough
%       correction on. This index corresponds to the channel directories
%       location in the cell array movieData.channels). If not
%       input, the user will be asked to select from the movie's channels.       
%   
%       ('BleedChannelIndex'-> Positive integer scalar or vector)
%       The integer indices of the channel(s) whose fluorophore is bleeding
%       through into the channel to be corrected. You will need to have a
%       bleedthrough coefficient for each of these channels.
%       If not input, user will be asked to select, unless batch mode is
%       enabled, in which case an error will be generated.
%
%       ('BleedCoefficients' -> Positive scalar or vector)
%       The bleedthrough coefficients for each bleed channel specified by
%       BleedImageChannels, in the same order. This should be the average
%       coefficient produced by calculateMovieBleedthrough.m
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed.
%
%
%
% Output:
%
%   movieData - the updated movieData object with the correction
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The corrected images are written to the directory specified by the
%   parameter OuptuDirectory, with each channel in a separate
%   sub-directory. They will be stored as bit-packed .tif files. 
%
% 
% Hunter Elliott
% 11/2009
% Revamped 6/2010
%

%% ------ Parameters ------- %%

pString = 'btc_'; %The string to prepend before the bleedthrough-corrected image directory & channel name
dName = 'bleedthrough_corrected_images_for_channel_';%String for naming the directories for each corrected channel

%% ----------- Input ------------ %%

%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous bleedthrough correction processes from this
%function
iProc = find(cellfun(@(x)(isa(x,'BleedthroughCorrectionProcess')),movieData.processes_),1);                          

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(BleedthroughCorrectionProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%% Initialization

nChan = numel(movieData.channels_);


%Ask the user for the image channel to correct if not input
if isempty(p.ChannelIndex)  
    if p.BatchMode
        error('In batch mode, you must specify the channel to perform bleedthrough correction on!')
    else
        %As the user to select a channel.
        p.ChannelIndex = selectMovieChannels(movieData,0,'Select the channel to bleedthrough correct:');        
    end    
elseif length(p.ChannelIndex) > 1
    error('Only one channel may be bleedthrough corrected at a time!')
end

if p.ChannelIndex > nChan || p.ChannelIndex < 1 || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel number specified! Check ChannelIndex input!!')
end


%Ask the user for the channels which bleed into the images to be corrected, if not input
if isempty(p.BleedChannelIndex)
    if p.BatchMode
        error('In batch mode, you must specify the channels containing the bleedthrough images!')
    else
        %As the user to select channel(s)        
        p.BleedChannelIndex = selectMovieChannels(movieData,1,...
            ['Select the channels which are bleeding through into channel ' num2str(p.ChannelIndex) ':']);        
    end
end

% If using the output of an existing process
if ~isempty(p.ProcessIndex)
    assert(isa(movieData.processes_{p.ProcessIndex},'ImageProcessingProcess'));
    
    % Check which channels have been processed by the process
    status= movieData.processes_{p.ProcessIndex}.checkChannelOutput;
    assert(all(status([p.ChannelIndex p.BleedChannelIndex])),...
        'The channel to be corrected, and the bleedthrough channels must all have been processed prior to bleedthrough correction!');
end


nBleed = length(p.BleedChannelIndex);

%Check/get the bleed coefficients
if isempty(p.BleedCoefficients)
    
    if p.BatchMode
        error('In batch mode, you must specify bleedthrough coefficients!')
    else
        %Ask the user to input the coefficients
        for j = 1:nBleed
            userInput = inputdlg(['Channel #' num2str(p.BleedChannelIndex(j)) ' bleedthrough coefficient:'],...
                                  'Bleedthrough Coefficient Input',1);
            if isempty(userInput)
                error('You must input a bleedthrough coefficient for each bleedthrough channel to continue!')                                
            end

            userInput = str2double(userInput);

            if ~isnumeric(userInput) || userInput <= 0 || ~isfinite(userInput) || ~isreal(userInput) %Yea, I'd like to see you try to get past all those checks!!
                error('Invalid bleedthrough coefficient! Must be a finite positive number!')            
            end

            p.BleedCoefficients(j) = userInput;

        end
    end
elseif length(p.BleedChannelIndex) ~= length(p.BleedCoefficients)
    error('You must specify the same number of bleedthrough coefficients and bleedthrough images!')
end

%Log the selected channels in the process
movieData.processes_{iProc}.setPara(p);


%% Bleedthrough correction
disp('Starting bleedthrough correction...')


%Retrieve the paths and names of the input images   
if isempty(p.ProcessIndex)
    imPaths = movieData.getChannelPaths;
    imNames = movieData.getImageFileNames;
else
    imPaths = movieData.processes_{p.ProcessIndex}.outFilePaths_;
    imNames = movieData.processes_{p.ProcessIndex}.getOutImageFileNames;
end

inDir = imPaths{1,p.ChannelIndex}; %input images to be corrected
inNames = imNames(1,p.ChannelIndex);
bleedImDir = imPaths(1,p.BleedChannelIndex);
bleedImNames = imNames(1,p.BleedChannelIndex);

%Log the input and bleed images in the process
movieData.processes_{iProc}.setInImagePath(p.ChannelIndex,inDir);
movieData.processes_{iProc}.setCorrectionImagePath(p.BleedChannelIndex,bleedImDir);

% Set the output directory
outDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex)]; %Corrected images
mkClrDir(outDir);
movieData.processes_{iProc}.setOutImagePath(p.ChannelIndex,outDir);


%% -------------- Apply bleedthrough correction ------------%%
%Applies the bleedthrough correction from above to each selected channel

disp('Applying bleedthrough correction to images...')

%Go through each image and apply the appropriate bleedthrough correction

if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, bleedthrough correcting channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end        

nImages = movieData.nFrames_;

disp(['Bleedthrough correcting channel ' num2str(p.ChannelIndex) '...']);
disp(['Correcting images from ' inDir]);
disp(['Storing results in ' outDir]);
for j = 1:nBleed
    disp(['Correcting bleedthrough from channel ' num2str(p.BleedChannelIndex(j)) ', using images from ' bleedImDir{j}])        
    disp(['...using bleedthrough coefficient of ' num2str(p.BleedCoefficients(j))])
end
for iImage = 1:nImages
    
    
    %Load the image to be corrected
    currIm = imread([inDir filesep inNames{1}{iImage}]);
    
    %Check the bit-depth of the image
    ogClass = class(currIm);
    currIm = double(currIm);

    for iBleed = 1:nBleed        
        
        %Load the bleed image
        currBleedIm = double(imread([bleedImDir{iBleed} filesep bleedImNames{iBleed}{iImage}]));

        %Subtract the bleedthrough from this channel
        currIm = currIm - (currBleedIm * p.BleedCoefficients(iBleed));
        
        %Remove negative values (these usually occur in the background)
        currIm(currIm < 0) = 0;
        
        if ~any(currIm > 0)
            error('Please check the specified bleedthrough coefficients: The specified correction results in completely blank images!')
        end

    end                

    %Cast to original class
    currIm = cast(currIm,ogClass);

    %Write it to disk    
    imwrite(currIm,[outDir filesep pString inNames{1}{iImage}]);

    if ~p.BatchMode && mod(iImage,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iImage / nImages,wtBar)
    end                        


end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------------- Output ------- %%

disp('Saving results...')

%Log the correction in the movieData object and save it

movieData.processes_{iProc}.setDateTime;
movieData.save;


disp('Finished Correcting!')