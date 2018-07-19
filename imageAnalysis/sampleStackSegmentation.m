function samples = sampleStackSegmentation(mask,stack,ind)
%SAMPLESTACKSEGMENTATION samples each input image in the stack in the segmented region(s)
%
% samples = sampleStackSegmentation(mask,stack)
% samples = sampleStackSegmentation(mask,stack,ind)
%
%   mask - MxN mask or index image specifying areas to segment. If logical,
%   only a single sample will be returned (run bwlabel first to get index
%   image)
%
%   stack - MxNxP stack of images to sample, where P is the number of
%   images.
%
%   ind - indices to sample. Optional. If not input, all unique non-zero
%   values will be sampled.
%

%Hunter Elliott
%4/2013

%% ------  Input  --------- %%

ip=inputParser;
ip.addRequired('mask',@ismatrix);
ip.addRequired('stack');
ip.parse(mask,stack);

assert(size(stack,1)==size(mask,1) && size(stack,2) == size(mask,2),...
    'The image and mask must have the same size along the first 2 dimensions!');


%% ----------- Init --------- %%

%Avoid rounding error on integer class inputs
stack = double(stack);

if nargin < 3 || isempty(ind)
    %Get indices
    ind = unique(mask(mask(:)>0));
end
nInd = numel(ind);

nIm = size(stack,3);


%Initialize sample array
sampledFields = {'avg','med','std','min','max','n'};%Field names
sampleFuns = {@mean,@median,@std,@min,@max,@numel};%Functions for calculating each field value
nFields = numel(sampledFields);
samples(size(stack,3),1)=struct();
for j = 1:nFields
    [samples.(sampledFields{j})] = deal(nan(nInd,1));
end

%Separate out the stack so we can do logical indexing below
stack = arrayfun(@(x)(stack(:,:,x)),1:nIm,'Unif',false);


%% ------------- Sampling ---------------- %%

for iInd = 1:nInd
   
    for iIm = 1:nIm
        
        currSamp = stack{iIm}(mask(:) == ind(iInd));
        
        for iField = 1:nFields                                                
            samples(iIm).(sampledFields{iField})(iInd) = sampleFuns{iField}(currSamp);    
        end
    end
    
end



