function [imfOut,noise] = testImf(imf,alpha)
%This function tests the IMF significance by comparing all imf points with
%the imf of a gaussian white noise
%It uses the Anderson-Darling test (thanks to Francois)
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
%Marco Vilela, 2012

[nImf,nObs] = size(imf);

%This numbers work well. They are just the size of the white noise
%genereted to test the imf's
nSurr    = 10;
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
        Ho = cell2mat(arrayfun(@(x) adtestk({imfTest,x},alpha),imf(i,:),'UniformOutput',0));
        imfOut(i,Ho == 1) = 0;
        noise   = noise + imfOut(i,:);
    end
end