function [imfOut,noise] = testImf(imf,varargin)
%This function tests the IMF significance by comparing all imf points with
%the imf of a gaussian white noise
%
%Usage:
%     [imfOut,noise] = testIMF(imf)
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
%       noise - estimated noise level
%
% Required functions: emdc, emd (extern)
%
%Reference:
%           The Hilbert-Huang Transform and Its Applications. Norden Huang,
%           Samuel Shen. World Scientific. Chapter 5.
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('imf',@ismatrix);
ip.addParamValue('alpha',0.05,@isscalar);
ip.addParamValue('nSurr',100,@isscalar);

ip.parse(imf,varargin{:});
alpha  = ip.Results.alpha;
nSurr  = ip.Results.nSurr;

[nImf,nObs] = size(imf);


%crude estimative of the noise variance using the first IMF (See reference)
surrP  = nObs;
wnStd  = std(imf(1,:));
Wn     = wnStd*randn(surrP,nSurr);%WN generated for the test
noise  = zeros(1,numel(imf(1,:)));
imfOut = imf;

%Decomposing noise into IMF's
flag_mexError = false;
for i=1:nSurr
    
    try 
        wnImf{i}  = emdc([],Wn(:,i)); %(T,X)
    catch
        flag_mexError = true;
%         warning('mex emdc failed, reverting to slower matlab-based emd');
        wnImf{i}  = emd(Wn(:,i), 'T', []);
    end
    
end

% if flag_mexError
%     disp('mex emdc failed, reverting to slower matlab-based emd');
% end

%Testing each imf
for i = 1:nImf
    
    imfTest = cell2mat(cellfun(@(x) getThImf(x,i),wnImf,'UniformOutput',0));
    %Test the ith imf against the noise ith imf
    if ~isempty(imfTest)
        
        limit  = prctile(imfTest,100*[alpha/2 (1 - alpha/2)]);
        %Null hypothesis that imf(i,j) is at the gaussian noise imf level
        Ho                = bsxfun(@gt,imf(i,:),limit(2)) | bsxfun(@lt,imf(i,:),limit(1));
        imfOut(i,Ho == 1) = 0;
        noise             = noise + imfOut(i,:);
        
    end
    
end

noise = noise(:);