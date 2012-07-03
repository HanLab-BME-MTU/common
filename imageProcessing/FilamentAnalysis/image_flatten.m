function movieData = image_flatten(movieData, varargin)

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('funParams',[], @isstruct);
ip.parse(movieData,varargin{:});
funParams=ip.Results.funParams;

selected_channels = funParams.ChannelIndex;
flatten_method_ind = funParams.method_ind;
Gaussian_sigma = funParams.GaussFilterSigma;

nFrame = movieData.nFrames_;
nChannel = length(MD.channels_);

ImageFlattenOutputDir = [movieData.outputDirectory_,'/ImageFlatten'];
if (~exist(ImageFlattenOutputDir,'dir'))
    mkdir(ImageFlattenOutputDir);
end

for iChannel = selected_channels
    
    % Make output directory for the flattened images
    ImageFlattenOutputDir = [movieData.outputDirectory_,'/ImageFlatten/Channel',num2str(iChannel)];
    if (~exist(ImageFlattenOutputDir,'dir'))
        mkdir(ImageFlattenOutputDir);
    end
    
    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImg = double(currentImg);
        
        % Normalize before flatten
        currentImg = currentImg - min(min(currentImg)) +1;
        currentImg = currentImg/(max(max(currentImg)));
        
        % based on the given method index, do log or sqrt to flatten the image
        if flatten_method_ind == 1
            currentImg = log(currentImg);
            currentImg = currentImg - min(min(currentImg)) +1;
            currentImg = currentImg/(max(max(currentImg)));
            
        else
            if flatten_method_ind == 2
                currentImg = sqrt(currentImg);
            end
        end

        % Smooth the image in requested
        if Gaussian_sigma > 0
            currentImg = imfilter(currentImg, fspecial('gaussian',round(5*Gaussian_sigma), Gaussian_sigma),'replicate','same');
        end
        imwrite(currentImg,[ImageFlattenOutputDir,'/channel',num2str(iChannel),'_flatten_',num2str(iFrame),'.tif']);
    end
end
