function [imfOut,noise] = testImf(imf,alpha)
%This function tests the IMF significance by comparing all imf points with
%the imf of a gaussian white noise
%
%Usage:
%     [auxImf,noise] = testIMF(imf)
%
%Input:
%      imf   - matrix with all imf's 
%      
%      alpha - significance level        
%
%Output:
%       imfOut - matrix same size as the input imf 
%                the signal is removed from imf -> set to zero
%
%       noise - estimated noise
%
% See also: emdc
%
%Reference:
%           The Hilbert-Huang Transform and Its Applications. Norden Huang,
%           Samuel Shen. World Scientific. Chapter 5.
%Marco Vilela, 2012

[nImf,nObs] = size(imf);

%This numbers work well. They are just the size of the white noise
%genereted to test the imf's
nSurr    = 100;
surrP    = nObs;

%crude estimative of the noise variance
wnStd  = std(imf(1,:));
Wn     = wnStd*randn(surrP,nSurr);%WN generated for the test
noise  = zeros(1,numel(imf(1,:)));
imfOut = imf;

%Decomposing noise into IMF's
for i=1:nSurr
    wnImf{i}  = emdc([],Wn(:,i));
end

%Testing each imf
for i = 1:nImf
    imfTest = cell2mat(cellfun(@(x) getThImf(x,i),wnImf,'UniformOutput',0));
    %Test the ith imf against the noise ith imf
    if ~isempty(imfTest)
        limit  = prctile(imfTest,100*[alpha/2 (1 - alpha/2)]);
        %Null hypothesis that imf(i,j) is at the gaussian noise imf level
        Ho     = bsxfun(@gt,imf(i,:),limit(2)) | bsxfun(@lt,imf(i,:),limit(1));
        imfOut(i,Ho == 1) = 0;
        noise  = noise + imfOut(i,:);
    end
end